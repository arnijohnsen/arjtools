#' Location of probes in hg18 385k design array
#'
#' A data.table with genomic coordinates of 385815 probes in an array,
#'   using hg18 assembly. Extracted from design files.
#'
#' @format A data table with 385815 obs. of 2 variables:
#' \describe{
#'   \item{chrom}{chromosome, in the format "chr1", "chrY", etc.}
#'   \item{pos}{position, position of probe on chromosome}
#' }
#' @source Desgin files
"hg18_385k_pos"

#' Location of probes in hg17 385k design array
#'
#' A data.table with genomic coordinates of 388560 probes in an array,
#'   using hg17 assembly. Extracted from design files.
#'
#' @format A data table with 388560 obs. of 2 variables:
#' \describe{
#'   \item{chrom}{chromosome, in the format "chr1", "chrY", etc.}
#'   \item{pos}{position, position of probe on chromosome}
#' }
#' @source Desgin files
"hg17_385k_pos"

#' Map from probes in hg17 385k array design to locations in hg18 assembly
#'
#' A data.table with genomic coordinates of 388560 probes in an array,
#'   using hg17 assembly. Each probe has listed the corresponding genomic
#'   coordinates in hg18 assembly.
#'
#' @format A data table with 388560 obs. of 4 variables:
#' \describe{
#'   \item{hg17chrom}{chromosome of probe in hg17, in the format "chr1", "chrY", etc.}
#'   \item{hg17pos}{position of probe in hg17, position of probe on chromosome}
#'   \item{hg18chrom}{chromosome of corresponding location in hg17, in the format "chr1", "chrY", etc.}
#'   \item{hg18pos}{position of corresponding location in hg17, on chromosome}
#' }
#' @source \url{https://genome.ucsc.edu/cgi-bin/hgLiftOver}
"hg17_hg18_map"
