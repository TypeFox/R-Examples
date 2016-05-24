context("Meaningful error on bad input data")

test_that("Number of observations in phenotype and snp matrix are equal",{
  famFile <- system.file("extdata", "plinkPhenotypeExample.fam", package = "geneSLOPE")
  mapFile <- system.file("extdata", "plinkMapExample.map", package = "geneSLOPE")
  snpsFile <- system.file("extdata", "plinkDataExample.raw", package = "geneSLOPE")
  phenotype <- read_phenotype(filename = famFile)

  phenotype$y <- rnorm(91)

  expect_error(screen_snps(snpsFile, mapFile, phenotype, pValMax = 0.05,
                        chunkSize = 1e2, verbose=FALSE), "phenotype")
})

test_that("p-value has to be smaller than 1, larger than 0",{
  famFile <- system.file("extdata", "plinkPhenotypeExample.fam", package = "geneSLOPE")
  mapFile <- system.file("extdata", "plinkMapExample.map", package = "geneSLOPE")
  snpsFile <- system.file("extdata", "plinkDataExample.raw", package = "geneSLOPE")
  phenotype <- read_phenotype(filename = famFile)

  expect_error(screen_snps(snpsFile, mapFile, phenotype, pValMax = 1.2,
                        chunkSize = 1e2, verbose=FALSE), "pValMax")
  expect_error(screen_snps(snpsFile, mapFile, phenotype, pValMax = -0.5,
                        chunkSize = 1e2, verbose=FALSE), "pValMax")
})

test_that("chunk size has to be positive",{
  famFile <- system.file("extdata", "plinkPhenotypeExample.fam", package = "geneSLOPE")
  mapFile <- system.file("extdata", "plinkMapExample.map", package = "geneSLOPE")
  snpsFile <- system.file("extdata", "plinkDataExample.raw", package = "geneSLOPE")
  phenotype <- read_phenotype(filename = famFile)

  expect_error(screen_snps(snpsFile, mapFile, phenotype, pValMax = 0.05,
                        chunkSize = -100, verbose=FALSE), "chunkSize.*positive")
})


test_that("Give warning when chunk size is smaller than 10",{
  famFile <- system.file("extdata", "plinkPhenotypeExample.fam", package = "geneSLOPE")
  mapFile <- system.file("extdata", "plinkMapExample.map", package = "geneSLOPE")
  snpsFile <- system.file("extdata", "plinkDataExample.raw", package = "geneSLOPE")
  phenotype <- read_phenotype(filename = famFile)

  expect_warning(screen_snps(snpsFile, mapFile, phenotype, pValMax = 0.05,
                        chunkSize = 9, verbose=FALSE), "chunk")
})

#
# test_that("Warning when parameter rawFile is not .raw",{
#   famFile <- system.file("extdata", "plinkPhenotypeExample.fam", package = "geneSLOPE")
#   mapFile <- system.file("extdata", "plinkMapExample.map", package = "geneSLOPE")
#   # we change file extension
#   snpsFile <- "plinkDataExample.rwa"
#   phenotype <- read_phenotype(filename = famFile)
#
#   expect_warning(screen_snps(snpsFile, mapFile, phenotype, pValMax = 0.05,
#                           chunkSize = 100, verbose=FALSE), ".raw")
# })


test_that("Meaningful error when file cannot be found",{
  famFile <- system.file("extdata", "plinkPhenotypeExample.fam", package = "geneSLOPE")
  mapFile <- system.file("extdata", "plinkMapExample.map", package = "geneSLOPE")
  # we change file to nonexisting
  snpsFile <- "nonExistingFile.raw"
  phenotype <- read_phenotype(filename = famFile)

  expect_error(screen_snps(snpsFile, mapFile, phenotype, pValMax = 0.05,
                          chunkSize = 100, verbose=FALSE), "not found")
})
