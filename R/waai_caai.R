#' Compute whole arm abberation index (waai) for a chromosome arm
#'
#' @param seg_values Numeric vector, with values associated with segments from one chromosome arm (usually mean log ratio)
#' @param seg_nprobes Numeric vector, with number of probes in each segment
#' @return Whole arm abberation index for the given chromosome arm
#' @examples
#' # Read aCGH data to data.frame dat
#' seg <- copynumber::pcf(dat)
#' seg.dt <- data.table::data.table(seg)
#' seg.dt[,.(waai = waai(mean, n.probes)), by=.(chrom, arm)]
#' @export
waai <- function(seg_values, seg_nprobes){
  return(weighted.mean(seg_values, seg_nprobes, na.rm = T))
}
