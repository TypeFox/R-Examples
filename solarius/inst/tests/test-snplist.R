
context("snplist/snpind")

# input data
test_that("input data parameters `genocovdata`", {
  data(dat50)
  snplist <- c("s1", "s2")
  
  mod <- solarAssoc(trait ~ 1, phenodata, snpdata = genodata, snplist = snplist)
  
  expect_equal(nrow(mod$snpf), length(snplist))
})

test_that("1 `genocov.files` & subset of SNPs in 1 `snplists.files`", {
  data(dat50)
  
  dir <- system.file("inst/extdata/solarAssoc/", package = "solarius")
  genocov.files <- file.path(dir, "snp.genocov")
  snplists.files <- file.path(dir, "snp.geno-list2")

  snplists.length <- length(readLines(snplists.files))
  
  snpind <- 1
  mod <- solarAssoc(trait ~ 1, phenodata, genocov.files = genocov.files, snplists.files = snplists.files, snpind = snpind)
  
  expect_equal(nrow(mod$snpf), length(snpind))
})

test_that("(error) `snpind` for 2 `genocov.files` & subset of SNPs in 2 `snplists.files`", {
  data(dat50)
  
  dir <- system.file("inst/extdata/solarAssoc/", package = "solarius")
  genocov.files <- file.path(dir, c("snp.genocov2", "snp.genocov3"))
  snplists.files <- file.path(dir, c("snp.geno-list2", "snp.geno-list3"))

  snpind <- 1:2  
  expect_error(mod <- solarAssoc(trait ~ 1, phenodata, genocov.files = genocov.files, snplists.files = snplists.files, snpind = snpind))
})
