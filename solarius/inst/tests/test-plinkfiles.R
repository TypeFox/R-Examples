
context("plinkfiles")

test_that("input `plink.raw`", {
  data(dat50)
  
  dir <- package.file("inst/extdata/solarAssoc/plink/", package = "solarius")
  plink.raw <- file.path(dir, "dat50.raw")
  
  num.snps <- 50
   
  mod <- solarAssoc(trait ~ 1, phenodata, plink.raw = plink.raw)
  expect_true(nrow(mod$snpf) == num.snps)
})

test_that("input `plink.raw` and `plink.map`", {
  data(dat50)
  
  dir <- package.file("inst/extdata/solarAssoc/plink/", package = "solarius")
  plink.raw <- file.path(dir, "dat50.raw")
  plink.map <- file.path(dir, "dat50.map")
  
  num.snps <- 50
   
  mod <- solarAssoc(trait ~ 1, phenodata, plink.raw = plink.raw, plink.map = plink.map)
  
  expect_true(nrow(mod$snpf) == num.snps)
  expect_true("pos" %in% colnames(mod$snpf))
  expect_true(sum(is.na(mod$snpf$pos)) == 0)
})

test_that("input `plink.ped` and `plink.map`", {
  data(dat50)
  
  mod0 <- solarAssoc(trait ~ 1, phenodata, snpdata = genodata)
  setkey(mod0$snpf, pSNP)
  snp.best0 <- mod0$snpf[1, SNP]

  # plink format  
  dir <- package.file("inst/extdata/solarAssoc/plink/", package = "solarius")
  plink.ped <- file.path(dir, "dat50.ped")
  plink.map <- file.path(dir, "dat50.map")
  
  num.snps <- 50
   
  mod <- solarAssoc(trait ~ 1, phenodata, plink.ped = plink.ped, plink.map = plink.map)
  setkey(mod$snpf, pSNP)
  snp.best <- mod$snpf[1, SNP]

  expect_true(nrow(mod$snpf) == num.snps)
  expect_true("pos" %in% colnames(mod$snpf))
  expect_true(sum(is.na(mod$snpf$pos)) == 0)

  expect_equal(snp.best0, snp.best)
})

