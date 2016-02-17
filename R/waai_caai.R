#' Compute whole-arm abberation index (waai) for a chromosome arm
#'
#' @param seg_values Numeric vector, with values associated with segments from one chromosome arm (usually mean log ratio).
#' @param seg_nprobes Numeric vector, with number of probes in each segment.
#' @return Whole-arm abberation index for the given chromosome arm.
#' @examples
#' # Read aCGH data to data.frame dat
#' seg <- copynumber::pcf(dat)
#' seg.dt <- data.table::data.table(seg)
#' seg.dt[,.(waai = waai(mean, n.probes)), by=.(chrom, arm)]
#' @export
waai <- function(seg_values, seg_nprobes){
  return(weighted.mean(seg_values, seg_nprobes, na.rm = T))
}

#' Compute complex arm-wise aberration index (caai) for a chromosome arm
#'
#' @param seg_values Numeric vector, with values associated with segments from one chromosome arm (usually mean log ratio).
#' @param seg_start Numeric vector, with genomic position of segment start.
#' @param seg_end Numeric vector, with genomic position of segment end.
#' @param alpha A number, alpha value for computation of P. Defaults to the value used by Russnes et al (2010).
#' @param thetaH A number, theta value for compuation of Q. Defaults to the value used by Russnes et al (2010).
#' @param R A number, size of window R over which to sum S. Defaults to the value used by Russnes et al (2010).
#' @return Complex arm-wise aberration index for the given chromosome arm.
#' @examples
#' # Read aCGH data to data.frame dat
#' seg <- copynumber::pcf(dat)
#' seg.dt <- data.table::data.table(seg)
#' seg.dt[,.(caai = caai(mean, start.pos, end.pos)), by=.(chrom, arm)]
#' @export
caai <- function(seg_values, seg_start, seg_end, alpha = 10000/0.005, thetaH = 1.2, R = 20e6){
  n <- length(seg_values)
  seg_lengths <- seg_end - seg_start + 1
  L1 <- seg_lengths[1:(n-1)]
  L2 <- seg_lengths[2:n]
  H1 <- seg_values[1:(n-1)]
  H2 <- seg_values[2:n]
  P <- tanh( alpha / ( L1 + L2 ))
  Q <- tanh( abs( H2 - H1 ) / thetaH )
  W <- 0.5*(1 + tanh( 10*(P - 0.5) / tanh(5) ) )
  S <- W*pmin(P, Q)
  sum_n <- findInterval(seg_start+R, seg_start)[1:(n-1)] - 1:(n-1) +1
  SR <- zoo::rollapply(S, sum_n, sum, align = "left", fill = NA)
  return(max(SR, na.rm = T))
}
