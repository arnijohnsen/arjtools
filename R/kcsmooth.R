#' Smooth samples with KCsmart algorithm
#'
#' Wrapper to smooth aCGH data with algorithm from KCsmart.
#'
#' @param dt data.table with aCGH data. First two colums should be names chrom and pos
#'   and contain information about the location of input probes. Chromosomes should be
#'   names "1", "Y", etc.
#' @param mirrorLocs List containing the chromosome start, centromere and end positions
#' @param sigma The kernel width
#' @param sampleDensity The sample point matrix resolution
#' @param verbose If set to false, no progress information is displayed
#' @param maxmem This parameter controls memory usage, set to lower value to lower memory consumption
#' @param what character, determines what data should be returned. "pos" for only positive KCscore,
#'   "neg" for only negative KCscore, "both" for both positive and negative KCscore.
#' @return data.table with smoothed values. First column is chromosome, second column indicates if value is from positive
#'   or negative smoothed values. All other columns are smoothed values for each sample.
#' @import data.table
#' @export
kcsmooth <- function(dt, mirrorLocs, sigma = 1e+06, sampleDensity = 50000, maxmem = 1000, verbose = T, what = "both"){
  if(!(what %in% c("pos", "neg", "both"))){
    stop("what must be pos, neg or both")
  }
  spm_to_dt <- function(spm, what){
    spm_unlist <- unlist(spm@data)
    dt <- data.table(
      chrom = str_replace(names(spm_unlist), "\\..*", ""),
      posneg  = str_extract(names(spm_unlist), "[a-z]+"),
      value = spm_unlist)
    if(what != "both"){
      dt <- dt[posneg == what]
    }
    return(dt)
  }
  # kcsmooth single data
  spm <- KCsmart::calcSpm(dt[,c(1:2, 3), with = F], hsMirrorLocs, sigma, sampleDensity, verbose)
  all_dt <- spm_to_dt(spm, what)
  n <- dim(dt)[2]
  if (n > 3){
    for(i in 4:n){
      spm <- KCsmart::calcSpm(dt[,c(1:2, i), with = F], hsMirrorLocs, sigma, sampleDensity, verbose)
      all_dt[,i:=spm_to_dt(spm, what)$value]
    }
  }
  setnames(all_dt, c("chrom", "what", names(dt)[3:n]))
  return(all_dt)
}
