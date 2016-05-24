
context("solarAssoc")

# input data
test_that("input data parameters `genocovdata` and `genodata`", {
  data(dat50)
  snps <- 1:2
  
  mod1 <- solarAssoc(trait ~ 1, phenodata, snpdata = genodata[, snps], kinship = kin)
  mod2 <- solarAssoc(trait ~ 1, phenodata, snpcovdata = genocovdata[, snps], kinship = kin)
  
  expect_equal(mod1$snpf$pSNP, mod2$snpf$pSNP)
})

test_that("input data parameters 2 `genocov.files` & 2 `snplists.files`", {
  data(dat50)
  
  dir <- system.file("inst/extdata/solarAssoc/", package = "solarius")
  genocov.files <- file.path(dir, c("snp.genocov2", "snp.genocov3"))
  snplists.files <- file.path(dir, c("snp.geno-list2", "snp.geno-list3"))

  mod <- solarAssoc(trait ~ 1, phenodata, genocov.files = genocov.files, snplists.files = snplists.files)
  expect_true("solarAssoc" %in% class(mod))
})
