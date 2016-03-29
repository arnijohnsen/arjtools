#' Finds nearest probe on hg18_385k
#'
#' Returns nearest probe on a 385k array design in hg18 assembly
#'
#' @param x Numeric vector, with genomic coordinate of probe
#' @param chrom Character or character vector (in the form "chr11", "chrX", etc.),
#'   with chromosome of probe. All elements must be identical.
#' @return Numeric vector, with genomic coordinate of closest probe in the hg18_385k design
#' @examples
#'   # Create sample data
#'   library(data.table)
#'   dat <- data.table(chrom = c("chr1", "chr1"), pos = c(46348, 46349))
#'   dat[, pos_385k := map_to_385k(pos, chrom), by = chrom]
#'   dat
#' @import data.table
#' @export
map_to_385k <- function(x, chrom){
  if(!all(chrom == chrom[1])){
    stop("All elements of chrom must be equal")
  }
  load(hg18_385k_pos)
  chr <- chrom[1]
  dt <- copy(hg18_385k_pos[chrom == chr])
  dt[,val := pos]
  setattr(dt, "sorted", "val")
  return(dt[J(x), roll = "nearest"]$pos)
}
