context("searchbyterm")

test_that("searchbyterm works correctly", {
  skip_on_cran()
  
  # Find multiple species
  a <- searchbyterm(genus = "ochotona", specificepithet = "(princeps OR collaris)", limit = 5, verbose = FALSE)
  b <- searchbyterm(class = "aves", state = "california", limit = 10, verbose = FALSE)
  cc <- searchbyterm(class = "aves", state = "california", year = 1976, limit = 10, verbose = FALSE)
  Sys.sleep(3)
  library("httr")
  #d <- searchbyterm(class = "aves", state = "california", year = ">1976", limit = 60, config = verbose())
  #d <- searchbyterm(class = "aves", state = "california", year = ">1976", limit = 60, verbose = FALSE)
  
  expect_equal( NROW( searchbyterm(limit = 1, verbose = FALSE)$data ), 1)
  expect_is(a, "list")
  expect_is(a$meta, "list")
  expect_is(a$data, "data.frame")
  expect_is(b$data, "data.frame")
  expect_is(cc, "list")
  #expect_is(d$meta, "list")
  
  expect_is(a$data$language, "character")
  #expect_match(d$data$class, "Aves")
  
  expect_equal(NROW(cc$data), 10)
  #expect_equal(NROW(d$data), 60)
  
  expect_equal(unique(as.numeric(cc$data$year)), 1976)
  
  #expect_gt(min(as.numeric(d$data$year)), 1976)
})


test_that("searchbyterm - state param works when using boolean's with > 1 state name", {
  skip_on_cran()
  
  aa <- searchbyterm(genus = "zapus", specificepithet = "hudsonius",  
                     stateprovince = "minnesota OR new mexico")
  
  expect_is(aa, "list")
  expect_is(aa$data, "data.frame")
  expect_equal(unique(tolower(aa$data$stateprovince)), "new mexico")
  
  aa <- searchbyterm(genus = "Ursus", stateprovince = "california OR florida", limit = 5)
  
  expect_is(aa, "list")
  expect_is(aa$data, "data.frame")
  expect_equal(unique(tolower(aa$data$stateprovince)), "california")
})


test_that("searchbyterm fails correctly", {
  skip_on_cran()
  
  # server error with inproper date
  expect_null(searchbyterm(specificepithet = "nigripes", date = "1935-09-abab", verbose = FALSE))
})

## FIX ME, cursor works when run alone, but not in this test
# test_that("searchbyterm cursor works correctly", {
#   out <- searchbyterm(genus = "Ochotona", limit = 1010)
#   expect_equal( NROW( out$data ), 1010)
#   expect_is(out, "list")
#   expect_is(out$meta, "list")
#   expect_is(out$data, "data.frame")
# })

test_that("searchbyterm multi-year param input works", {
  skip_on_cran()
  
  out <- suppressMessages(searchbyterm(gen = "ochotona", specificepithet = "(princeps OR collaris)", 
                      year = c(">=1916", "<=1920")))
  dates <- as.Date(na.omit(out$data$eventdate))
  asnumdates <- as.numeric(format(dates, "%Y"))
  expect_is(dates, "Date")
  expect_is(asnumdates, "numeric")
  expect_equal(min(asnumdates), 1916)
  expect_equal(max(asnumdates), 1920)
})
