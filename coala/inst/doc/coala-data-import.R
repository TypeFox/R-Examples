## ----import_data, results="hide"-----------------------------------------
suppressPackageStartupMessages(library(PopGenome))
fasta <- system.file("example_fasta_files", package = "coala")
data_pg <- readData(fasta, progress_bar_switch = FALSE)
data_pg <- set.outgroup(data_pg, c("Individual_Out-1", "Individual_Out-2"))
individuals <- list(paste0("Individual_1-", 1:5), paste0("Individual_2-", 1:5))
data_pg <- set.populations(data_pg, individuals)

## ------------------------------------------------------------------------
library(coala)
segsites <- as.segsites(data_pg)

## ----calc_sumstats_1-----------------------------------------------------
model <- coal_model(c(5, 5, 2), 1, 25) + 
  feat_mutation(5) +
  feat_outgroup(3) +
  sumstat_sfs(population = 1)
stats <- calc_sumstats_from_data(model, segsites)
stats$sfs

## ----calc_sumstats_2-----------------------------------------------------
stats <- calc_sumstats_from_data(model, data_pg)
stats$sfs

