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
probe_to_385k <- function(x, chrom){
  if(!all(chrom == chrom[1])){
    stop("All elements of chrom must be equal")
  }
  load(system.file("data", "hg18_385k_pos.rda", package = "arjtools"))
  chr <- chrom[1]
  dt <- copy(hg18_385k_pos[chrom == chr])
  dt[,val := pos]
  setattr(dt, "sorted", "val")
  return(dt[J(x), roll = "nearest"]$pos)
}

#' Converts aCGH data to hg18 385k design
#'
#' Converts aCGH data to a hg18 385k design.
#'   Each input probe is mapped to its nearest output probe in the 385k design.
#'   If multiple inputs match to a single output, their mean is computed. If no input
#'   matches to a 385k probe, an NA value is returned for that output.
#'
#' @param dt data.table with aCGH data. First two colums should be names chrom and pos
#'   and contain information about the location of input probes. Chromosomes should be
#'   names "chr1", "chrY", etc.
#'   The remaining columns are interpreted as samples with the column name denoting
#'   the sample name.
#' @param assembly String specifying assembly of input data. Can take values "hg18"
#'   (default) or "hg17" (in which case, design must match hg17 385k design).
#' @return data.table with mapped data. First two colums are chrom and pos
#'   and contain information about the location of input probes. The remaining columns
#'   are mapped data for each sample.
#' @import data.table
#' @export
dt_to_385k <- function(dt, assembly = "hg18"){
  if(!identical(names(dt)[1:2],c("chrom", "pos"))){
    stop("First two columns should be called chrom and pos")
  }
  if(!(assembly %in% c("hg17", "hg18"))){
    stop("Assembly must be either hg17 or hg18")
  }
  if(assembly == "hg17"){
    load(system.file("data", "hg17_hg18_map.rda", package = "arjtools"))
    dt$chrom <- hg17_hg18_map$hg18chrom
    dt$pos   <- hg17_hg18_map$hg18pos
    dt <- dt[!is.na(chrom)]
  }
  load(system.file("data", "hg18_385k_pos.rda", package = "arjtools"))
  dt[,pos := probe_to_385k(pos, chrom), by = chrom]
  setkey(dt, chrom, pos)
  return(copy(dt[,lapply(.SD, mean, na.rm = T), by = .(chrom, pos)][.(hg18_385k_pos$chrom, hg18_385k_pos$pos)]))
}
