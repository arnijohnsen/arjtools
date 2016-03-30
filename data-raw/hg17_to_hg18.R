library(data.table)

probes_hg17 <- fread(system.file("extdata", "probes_hg17_385k.txt", package = "arjtools"))
probes_hg18 <- fread(system.file("extdata", "probes_hg18_385k.txt", package = "arjtools"))
hg17_385k_pos <- fread(
  system.file("extdata", "hg17_385k_pos.txt", package = "arjtools"),
  select = c("CHROMOSOME", "POSITION")
)
setkey(hg17_385k_pos, CHROMOSOME, POSITION)
dt <- data.table(
  hg17chrom = probes_hg17$V1,
  hg17pos = probes_hg17$V2,
  hg18chrom = probes_hg18$V1,
  hg18pos = probes_hg18$V2
)
setkey(dt, hg17chrom, hg17pos)
hg17_hg18_map <- dt[.(hg17_385k_pos$CHROMOSOME, hg17_385k_pos$POSITION)]
devtools::use_data(hg17_hg18_map)
