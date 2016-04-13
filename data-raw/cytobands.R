library(data.table)

assemblies <- c("hg38", "hg19", "hg18", "hg17", "hg16")
cytoband_map <- list()
for(i in assemblies){
  dest <- system.file("extdata", str_c("cytoband_", i, ".txt.gz"), package = "arjtools")
  cytoband_map[[i]] <- fread(str_c("zcat ", dest))
}
devtools::use_data(cytoband_map)
