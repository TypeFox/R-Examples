
context("solarMultipoint")

# basic examples
test_that("basic example on dat30", {
  dat <- loadMulticPhen()
  mibddir <- package.file("extdata", "solarOutput", "solarMibds", package = "solarius")  
  
  chr <- 5
  num.files <- length(list.files(mibddir, paste0("mibd.", chr)))

  mod <- solarMultipoint(trait1 ~ 1, dat, mibddir = mibddir, chr = chr)
  
  expect_equal(num.files, nrow(mod$lodf))
  expect_true(all(mod$lodf$LOD > 10))
})

test_that("CSV IBD matices", {
  dat <- loadMulticPhen()
  mibddir <- package.file("extdata", "solarOutput", "solarMibdsCsv", package = "solarius")  

  chr <- 5
  num.files <- length(list.files(mibddir, paste0("mibd.", chr)))
  
  mod <- solarMultipoint(trait1 ~ 1, dat, mibddir = mibddir, chr = chr)
  
  expect_equal(num.files, nrow(mod$lodf))
  expect_true(all(mod$lodf$LOD > 10))
})

# multipoint options
test_that("interval", {
  dat <- loadMulticPhen()
  mibddir <- package.file("extdata", "solarOutput", "solarMibdsCsv", package = "solarius")  

  mod <- solarMultipoint(trait1 ~ 1, dat, mibddir = mibddir, chr = 5, interval = 5, multipoint.settings = "finemap off")
  
  expect_equal(nrow(mod$lodf), 2)
})

# IDs in IBDs/phen
test_that("ids(phen) is a subset in ids(IBD)", {
  data(dat30)
  mibddir <- package.file("extdata", "solarOutput", "solarMibdsCsv", package = "solarius")  

  mod <- solarMultipoint(trait1 ~ 1, dat30, mibddir = mibddir, chr = 5, interval = 5, multipoint.settings = "finemap off")
  
  expect_equal(nrow(mod$lodf), 2)
})

test_that("ids(IBD) is a subset in ids(phen)", {
  data(dat30)
  mibddir <- package.file("extdata", "solarOutput", "solarMibdsCsvIncomplete", package = "solarius")  

  suppressWarnings({
    mod <- solarMultipoint(trait1 ~ 1, dat30, mibddir = mibddir)
  })

  expect_equal(nrow(mod$lodf), 2)  
  expect_true(all(mod$lodf$LOD > 1.6))
})

# Bivariate linkage
test_that("bivariate linkage", {
  data(dat30)
  mibddir <- package.file("extdata", "solarOutput", "solarMibdsCsv", package = "solarius")  

  mod <- solarMultipoint(trait1 + trait2 ~ 1, dat30, mibddir = mibddir, chr = 5, interval = 5, multipoint.settings = "finemap off")
  
  expect_equal(nrow(mod$lodf), 2)
  expect_true(all(mod$lodf$LOD > 0.9))
})

# Compute linkage in parallel
test_that("linkage in parallel", {
  data(dat30)
  mibddir <- package.file("extdata", "solarOutput", "solarMibdsCsv", package = "solarius")  

  mod <- solarMultipoint(trait1 ~ 1, dat30, mibddir = mibddir, interval = 5, chr = c(2, 5), multipoint.settings = "finemap off", cores = 2)
  
  expect_equal(nrow(mod$lodf), 4)
  expect_true(all(mod$lodf$LOD > 1.4))
})

