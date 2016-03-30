#' Genomic coordinate to chromosome arm
#'
#' Returns chromosome arms for given chromosome and genomic position.
#'   Currently not implemented and returns NULL.
#'
#' @param chromosome Character or numeric vector, with chromosome of genomic coordinate
#' @param position Numeric vector, with genomic position within chromosome
#' @param assembly a string specifying which genome assembly version should be applied
#'   to determine chromosome arms. Allowed options are "hg38", hg19", "hg18", "hg17"
#'   and "hg16" (corresponding to the five latest human genome annotations in the
#'   UCSC genome browser).
#' @return Character vector, with choromosome arm of given genomic coordinates
coord_to_arm <- function(chromosome, position, assembly = "hg19"){
  if (!(assembly %in% c("hg38", "hg19", "hg18", "hg17", "hg16"))) {
    stop("Invalid assembly Allowed options are hg38, hg19, hg18, hg17 and hg16")
  }
  return(NULL)
}
