# Environment dataDir ----------------------------------------------------------

testthat::context("Environment dataDir")

testthat::test_that("setSpatialDataDir and getSpatialDataDir work correctly", {
  setSpatialDataDir("~")
  testthat::expect_equal(path.expand("~"), getSpatialDataDir())
  setSpatialDataDir(getwd())
  testthat::expect_equal(getwd(), getSpatialDataDir())
})

# 'get' function test ----------------------------------------------------------

testthat::context("'get' function test")

testthat::test_that("get functions handle errors correctly", {
  
  testthat::expect_error(getCountryCode(), 
               'argument "lon" is missing, with no default')
  testthat::expect_warning(getCountryCode(0,0), 
                 "1 locations appear to be over international waters and no countryCode can be assigned")
  testthat::expect_warning(getCountryCode(c(0,0), c(0,45)), 
                 "1 locations appear to be over international waters and no countryCode can be assigned")
  testthat::expect_error(getCountryCode(0,100))
  
})

testthat::test_that("get functions return correct name", {
  
  testthat::expect_match(getCountryCode(2, 47), "FR")
  testthat::expect_match(getCountryCode(-80, 40), "US")  
  testthat::expect_match(getCountryCode(c(120,-17), c(-1.5,15)), "[IN SN]")
  
  testthat::expect_match(getTimezone(2, 47), "Europe/Paris")
  testthat::expect_match(getTimezone(-80, 40), "America/New_York")  
  testthat::expect_match(getTimezone(c(120,-17), c(-1.5,15)), "[Asia/Makassar Africa/Dakar]")
  
})

testthat::test_that("subsetting with countryCodes works", {
  
  testthat::expect_match(getCountryCode(2, 47), "FR")
  testthat::expect_match(getCountryCode(2, 47, countryCodes=c("FR")), "FR")
  testthat::expect_match(getCountryCode(2, 47, countryCodes="FR"), "FR")
  testthat::expect_warning(getCountryCode(2, 47, countryCodes=c("US")), 
                 "1 locations appear to be either over international waters or not in given countryCodes and no countryCode can be assigned")
  
  testthat::expect_match(getTimezone(2, 47), "Europe/Paris")
  testthat::expect_match(getTimezone(2, 47, countryCodes=c("FR")), "Europe/Paris")
  testthat::expect_match(getTimezone(2, 47, countryCodes="FR"), "Europe/Paris")
  testthat::expect_warning(getTimezone(2, 47, countryCodes=c("US")), 
                 "1 locations appear to be either over international waters or not in given countryCodes and no timezone can be assigned")
  
})

testthat::test_that("allData returns are correct dimension and type", {
  
  testthat::expect_is(getCountryCode(2, 47, allData=TRUE), "data.frame")
  testthat::expect_equal(dim(getCountryCode(2, 47, allData=TRUE)), c(1,11))
  testthat::expect_is(getCountryCode(c(120,-17), c(-1.5,15), allData=TRUE), "data.frame")
  testthat::expect_equal(dim(getCountryCode(c(120,-17), c(-1.5,15), allData=TRUE)), c(2,11))
  
  testthat::expect_is(getTimezone(2, 47, allData=TRUE), "data.frame")
  testthat::expect_equal(dim(getTimezone(2, 47, allData=TRUE)), c(1,6))
  testthat::expect_is(getTimezone(c(120,-17), c(-1.5,15), allData=TRUE), "data.frame")
  testthat::expect_equal(dim(getTimezone(c(120,-17), c(-1.5,15), allData=TRUE)), c(2,6))
  
})
