library(snpStats)

### read data
pedfile <- system.file("extdata/sample.ped.gz", package = "snpStats")
infofile <- system.file("extdata/sample.info", package = "snpStats")

sample <- read.pedfile(pedfile, snps = infofile)

### data tables
phen <- sample$fam

phen <- within(phen, {
  id <- rownames(phen)
  affected[is.na(affected)] <- 1
})
  
