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
#' @param groups Numeric vector. If NULL is passed (default), each sample is smoothed seperately.
#'   If a numeric vector is passed, it indicates which samples should be grouped together before
#'   smoothing. For example, groups = c(1,1,2) indicates that the first two samples will be grouped
#'   before smothing and the third will be smoothed seperately.
#'   If a non-NULL value is passed, the return data frame will have columns equal to the number of groups passed.
#'   (This is the same as passing group = c(1,2,3,...)).
#' @return data.table with smoothed values. First column is chromosome, second column indicates if value is from positive
#'   or negative smoothed values. All other columns are smoothed values for each sample or group.
#' @import data.table
#' @export
kcsmooth <- function(dt, mirrorLocs, sigma = 1e+06, sampleDensity = 50000, maxmem = 1000, verbose = T, what = "both", groups= NULL){
  if(!(what %in% c("pos", "neg", "both"))){
    stop("what must be pos, neg or both")
  }
  n <- dim(dt)[2]
  if(is.null(groups)){
    groups <- 1:(n-2)
  }
  n_groups <- length(unique(groups))
  spm_to_dt <- function(spm, what){
    spm_unlist <- unlist(spm@data)
    dt <- data.table(
      chrom = stringr::str_replace(names(spm_unlist), "\\..*", ""),
      posneg  = stringr::str_extract(names(spm_unlist), "[a-z]+"),
      value = spm_unlist)
    if(what != "both"){
      dt <- dt[posneg == what]
    }
    return(dt)
  }
  # Smooth first group
  spm <- KCsmart::calcSpm(dt[,c(1:2, which(groups == 1)+2), with = F], hsMirrorLocs, sigma, sampleDensity, verbose)
  all_dt <- spm_to_dt(spm, what)
  # Smooth groups 2, 3, ...
  if (n_groups > 1){
    for(i in 2:n_groups){
      spm <- KCsmart::calcSpm(dt[,c(1:2, which(groups == i)+2), with = F], hsMirrorLocs, sigma, sampleDensity, verbose)
      all_dt[,paste("col", i, sep = ""):=spm_to_dt(spm, what)$value]
    }
  }

  if(identical(groups, 1:(n-2))){
    setnames(all_dt, c("chrom", "what", names(dt)[3:n]))
  } else {
    setnames(all_dt, c("chrom", "what", paste("group", 1:n_groups, sep = "")))
  }
  return(all_dt)
}
