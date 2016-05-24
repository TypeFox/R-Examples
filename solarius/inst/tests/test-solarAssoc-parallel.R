
context("solarAssoc in parallel")

#----------------------------
# input data: `genocov.files` / `snpcovdata` / ...
#----------------------------
test_that("input `genocov.files`", {
  data(dat50)
  
  dir <- package.file("inst/extdata/solarAssoc/", package = "solarius")
  genocov.files <- file.path(dir, "snp.genocov")
  snplists.files <- file.path(dir, "snp.geno-list")

  mod1 <- solarAssoc(trait ~ age + sex, phenodata, genocov.files = genocov.files, snplists.files = snplists.files)
  mod2 <- solarAssoc(trait ~ age + sex, phenodata, genocov.files = genocov.files, snplists.files = snplists.files, cores = 2)
  
  mod3 <- solarAssoc(trait ~ 1, phenodata, genocov.files = genocov.files, snplists.files = snplists.files, cores = 2, batch.size = 11)
    
  expect_true(modelParCPUtime(mod1) > modelParCPUtime(mod2))
  expect_equal(modelParNumBatches(mod3), 5) # ceiling(50 SNPs / 11 batch size)
})

test_that("input `snpcovdata`", {
  data(dat50)
  
  mod <- solarAssoc(trait ~ 1, phenodata, snpcovdata = genocovdata, cores = 2, batch.size = 10)
  
  expect_equal(modelParNumBatches(mod), 5) # ceiling(50 SNPs / 11 batch size)
})

test_that("input `genodata`", {
  data(dat50)
  
  mod <- solarAssoc(trait ~ 1, phenodata, snpdata = genodata, cores = 2, batch.size = 10)
  
  expect_equal(modelParNumBatches(mod), 5) # ceiling(50 SNPs / 11 batch size)
})
