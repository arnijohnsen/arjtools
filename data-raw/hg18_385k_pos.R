library(data.table)

hg18_385k_pos <- fread(
  system.file("extdata", "hg18_385k_pos.txt", package = "arjtools"),
  select = c("CHROMOSOME", "POSITION")
)
setnames(hg18_385k_pos, c("chrom", "pos"))
setkey(hg18_385k_pos, chrom, pos)
devtools::use_data(hg18_385k_pos)

hg17_385k_pos <- fread(
  system.file("extdata", "hg17_385k_pos.txt", package = "arjtools"),
  select = c("CHROMOSOME", "POSITION")
)
setnames(hg17_385k_pos, c("chrom", "pos"))
setkey(hg17_385k_pos, chrom, pos)
devtools::use_data(hg17_385k_pos)
