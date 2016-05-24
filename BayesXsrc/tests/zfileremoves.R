## remove generated BayesX output files
testfiles <- c("mcmc.prg", "reml.prg", "step.prg",
  "mcmc.R", "reml.R", "step.R",
  "BayesX-tests.R", "data.raw")
files <- list.files()
files <- files[!files %in% testfiles]
files <- files[!grepl(".Rout", files)]
files <- files[!grepl(".save", files)]
file.remove(files)
