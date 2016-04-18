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
coord_to_arm <- function(chromosome, position, assembly = "hg19", full = F){
  if(length(chromosome) !=  length(position)){
    stop("chromosome and position must have equal length")
  }
  if (!(assembly %in% c("hg38", "hg19", "hg18", "hg17", "hg16"))) {
    stop("Invalid assembly, allowed options are hg38, hg19, hg18, hg17 and hg16")
  }
  if(stringr::str_sub(chromosome, 1, 3) != "chr"){
    chromosome <- str_c("chr", chromosome)
  }
  if(!grepl("chr[X-Y]|[0-9]+", chromosome)){
    stop("Invalid chromosome, must be 1-22, X or Y (or chr1-chr22, chrX or chrY)")
  }
  data(cytoband_map)
  names(position) <- chromosome
  arms <- rep("", length(chromosome))
  for(i in unique(chromosome)){
    map <- cytoband_map[[assembly]][V1 == i]
    arm <- map[(findInterval(position[i], map$V3)+1)]$V4
    if(!full){
      arm <- str_sub(arm, 1,1)
    }
    arms[chromosome == i] <- arm
  }
  return(arm)
}
