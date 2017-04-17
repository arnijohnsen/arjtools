#' Smooth samples with KCsmart algorithm
#'
#' Wrapper to smooth aCGH data with algorithm from KCsmart.
#'
#' @param dat data.table with aCGH data. First two colums should be names chrom and pos
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
kcsmooth <- function(dat, mirrorLocs, sigma = 1e+06, sampleDensity = 50000, maxmem = 1000, verbose = T, what = "both", groups = NULL){
  if(!(what %in% c("pos", "neg", "both"))){
    stop("what must be pos, neg or both")
  }
  if(!(all(names(dat)[1:2] == c("chrom", "maploc")))){
    setnames(dat, names(dat)[1:2], c("chrom", "maploc"))
  }
  if(stringr::str_detect(dat[1]$chrom, "^chr")){
    dat[,chrom := gsub("chr", "", chrom)]
  }
  n <- dim(dat)[2]
  if(is.null(groups)){
    groups <- 1:(n-2)
  }
  n_groups <- length(unique(groups))
  spm_to_dat <- function(spm, what){
    spm_unlist <- unlist(spm@data)
    dat <- data.table(
      chrom = stringr::str_replace(names(spm_unlist), "\\..*", ""),
      posneg  = stringr::str_extract(names(spm_unlist), "[a-z]+"),
      num = stringr::str_extract(names(spm_unlist), "[0-9]+$"),
      value = spm_unlist)
    if(what != "both"){
      dat <- dat[posneg == what]
    }
    return(dat)
  }
  # Smooth first group
  spm <- KCsmart::calcSpm(dat[,c(1:2, which(groups == 1)+2), with = F], hsMirrorLocs, sigma, sampleDensity, verbose)
  all_dat <- spm_to_dat(spm, what)
  # Smooth groups 2, 3, ...
  if (n_groups > 1){
    for(i in 2:n_groups){
      spm <- KCsmart::calcSpm(dat[,c(1:2, which(groups == i)+2), with = F], hsMirrorLocs, sigma, sampleDensity, verbose)
      all_dat[,paste("col", i, sep = ""):=spm_to_dat(spm, what)$value]
    }
  }
  if(identical(groups, 1:(n-2))){
    setnames(all_dat, c("chrom", "what", "num", names(dat)[3:n]))
  } else {
    setnames(all_dat, c("chrom", "what", "num", paste("group", 1:n_groups, sep = "")))
  }
  return(all_dat)
}

#' Create index from sig_region object
#'
#' Function to extract information from a sig_region object, to a list of all kcprobes within
#' those regions
#'
#' @param sig_regions A sig_regions object from kcsmart
#' @return data.table with 4 columns: number of gain/loss, chrom(osome) of region, whether its a gain or loss (pos/neg)
#'   and the num(ber) of kcprobe (genomic location of kcprobes is just sampleDensity*(number-1))
#' @import data.table
#' @export
kcmakefilter <- function(sig_regions){
  dat_list <- list()
  n_pos <- length(sig_regions@gains)
  for(i in 1:n_pos){
    dat_list[[i]] <- data.table(region = i,
                                chrom = sig_regions@gains[[i]]$chromosome,
                                what  = "pos",
                                num   = sig_regions@gains[[i]]$x)
  }
  n_neg <- length(sig_regions@losses)
  for(i in 1:n_neg){
    dat_list[[i+n_pos]] <- data.table(region = i,
                                chrom = sig_regions@losses[[i]]$chromosome,
                                what  = "neg",
                                num   = sig_regions@losses[[i]]$x)
  }
  dat <- rbindlist(dat_list)
  dat$num <- as.numeric(dat$num)
  return(dat)
}
