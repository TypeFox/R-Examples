
context("mgafiles")

test_that("input data parameters 2 `genocov.files` & 2 `snplists.files`", {
  data(dat50)
  
  dir <- system.file("inst/extdata/solarAssoc/", package = "solarius")
  genocov.files <- file.path(dir, c("snp.genocov2", "snp.genocov3"))
  snplists.files <- file.path(dir, c("snp.geno-list2", "snp.geno-list3"))
  
  snplists.length <- sum(laply(snplists.files, function(x) length(readLines(x))))
  
  mga.files <- list(genocov.files = genocov.files, snplists.files = snplists.files)
  
  mod <- solarAssoc(trait ~ 1, phenodata, mga.files = mga.files)
  expect_true(nrow(mod$snpf) == snplists.length)
})

test_that("`2 genocov.files` & 2 `snplists.files` & 2 `snpmap.files`", {
  data(dat50)
  
  dir <- system.file("inst/extdata/solarAssoc/", package = "solarius")
  genocov.files <- file.path(dir, c("snp.genocov2", "snp.genocov3"))
  snplists.files <- file.path(dir, c("snp.geno-list2", "snp.geno-list3"))
  snpmap.files <- file.path(dir, c("map.snp.1.2", "map.snp.1.3"))

  snplists.length <- sum(laply(snplists.files, function(x) length(readLines(x))))

  mga.files <- list(genocov.files = genocov.files, snplists.files = snplists.files, snpmap.files = snpmap.files)
  
  mod <- solarAssoc(trait ~ 1, phenodata, mga.files = mga.files)
  
  expect_true(nrow(mod$snpf) == snplists.length)
  expect_true("pos" %in% colnames(mod$snpf))
  expect_true(sum(is.na(mod$snpf$pos)) == 0)
})
