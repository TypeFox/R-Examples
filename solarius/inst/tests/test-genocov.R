
context("genocov")

test_that("1 `genocov.files` & subset of SNPs in 1 `snplists.files`", {
  data(dat50)
  
  dir <- system.file("inst/extdata/solarAssoc/", package = "solarius")
  genocov.files <- file.path(dir, "snp.genocov")
  snplists.files <- file.path(dir, "snp.geno-list2")

  snplists.length <- length(readLines(snplists.files))
  
  mod <- solarAssoc(trait ~ 1, phenodata, genocov.files = genocov.files, snplists.files = snplists.files)
  
  expect_equal(nrow(mod$snpf), snplists.length)
})


test_that("`1 genocov.files` & subset of SNPs in 2 `snplists.files`", {
  data(dat50)
  
  dir <- system.file("inst/extdata/solarAssoc/", package = "solarius")
  genocov.files <- file.path(dir, "snp.genocov")
  snplists.files <- file.path(dir, c("snp.geno-list2", "snp.geno-list3"))

  snplists.length <- sum(laply(snplists.files, function(x) length(readLines(x))))
  
  mod <- solarAssoc(trait ~ 1, phenodata, genocov.files = genocov.files, snplists.files = snplists.files)
  
  expect_equal(nrow(mod$snpf), snplists.length)
})

test_that("`1 genocov.files` & subset of SNPs in 2 `snplists.files` (2 cores)", {
  data(dat50)
  
  dir <- system.file("inst/extdata/solarAssoc/", package = "solarius")
  genocov.files <- file.path(dir, "snp.genocov")
  snplists.files <- file.path(dir, c("snp.geno-list2", "snp.geno-list3"))

  snplists.length <- sum(laply(snplists.files, function(x) length(readLines(x))))
  
  mod <- solarAssoc(trait ~ 1, phenodata, genocov.files = genocov.files, snplists.files = snplists.files, cores = 2)
  
  expect_equal(nrow(mod$snpf), snplists.length)
})


test_that("(error) `2 genocov.files` & 1 `snplists.files`", {
  data(dat50)
  
  dir <- system.file("inst/extdata/solarAssoc/", package = "solarius")
  genocov.files <- file.path(dir, c("snp.genocov2", "snp.genocov3"))
  snplists.files <- file.path(dir, "snp.geno-list2")

  expect_error(mod <- solarAssoc(trait ~ 1, phenodata, genocov.files = genocov.files, snplists.files = snplists.files))
})

test_that("`2 genocov.files` & subset of SNPs in 2 `snplists.files`", {
  data(dat50)
  
  dir <- system.file("inst/extdata/solarAssoc/", package = "solarius")
  genocov.files <- file.path(dir, c("snp.genocov2", "snp.genocov3"))
  snplists.files <- file.path(dir, c("snp.geno-list2", "snp.geno-list3"))

  snplists.length <- sum(laply(snplists.files, function(x) length(readLines(x))))
  
  mod <- solarAssoc(trait ~ 1, phenodata, genocov.files = genocov.files, snplists.files = snplists.files)
  
  expect_equal(nrow(mod$snpf), snplists.length)
})
