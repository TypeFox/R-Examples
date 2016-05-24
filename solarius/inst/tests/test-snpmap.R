
context("snpmap")

# input data
test_that("3 snps, 2 snps in map", {
  data(dat50)
  
  snpind <- 1:3
  snpind.map <- 1:2
  mod <- solarAssoc(trait ~ 1, phenodata, snpdata = genodata, snpind = snpind, snpmap = snpdata[snpind.map, ])
  
  expect_equal(nrow(mod$snpf), length(snpind))
  expect_true("pos" %in% colnames(mod$snpf))
  expect_true(sum(!is.na(mod$snpf$pos)) == length(snpind.map))
})

test_that("1 `genocov.files`", {
  data(dat50)
  
  dir <- system.file("inst/extdata/solarAssoc/", package = "solarius")
  genocov.files <- file.path(dir, "snp.genocov")
  snplists.files <- file.path(dir, "snp.geno-list2")
  snpmap.files <- file.path(dir, "map.snp.1")
  
  snplists.length <- length(readLines(snplists.files))
  
  mod <- solarAssoc(trait ~ 1, phenodata, genocov.files = genocov.files, snplists.files = snplists.files, snpmap.files = snpmap.files)
  
  expect_true(nrow(mod$snpf) == snplists.length)
  expect_true("pos" %in% colnames(mod$snpf))
  expect_true(sum(is.na(mod$snpf$pos)) == 0)
})

test_that("`2 genocov.files` & 2 `snplists.files` & 2 `snpmap.files`", {
  data(dat50)
  
  dir <- system.file("inst/extdata/solarAssoc/", package = "solarius")
  genocov.files <- file.path(dir, c("snp.genocov2", "snp.genocov3"))
  snplists.files <- file.path(dir, c("snp.geno-list2", "snp.geno-list3"))
  snpmap.files <- file.path(dir, c("map.snp.1.2", "map.snp.1.3"))
  
  snplists.length <- sum(laply(snplists.files, function(x) length(readLines(x))))
  
  mod <- solarAssoc(trait ~ 1, phenodata, genocov.files = genocov.files, snplists.files = snplists.files, snpmap.files = snpmap.files)
  
  expect_true(nrow(mod$snpf) == snplists.length)
  expect_true("pos" %in% colnames(mod$snpf))
  expect_true(sum(is.na(mod$snpf$pos)) == 0)
})
