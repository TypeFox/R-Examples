context("bigsearch")

test_that("bigsearch sends email correctly", {
  skip_on_cran()
  
  a <- bigsearch(specificepithet = "princeps", genus = "ochotona", 
                 rf = "test-bigsearch1", email = "myrmecocystus@gmail.com", verbose=FALSE)
  
	expect_null(a)
	expect_null(bigsearch(specificepithet = "collaris", genus = "ochotona", year = 1960,
	                      rf = "test-bigsearch2", email = "myrmecocystus@gmail.com", verbose=FALSE))
})

test_that("bigsearch fails correctly", {
  skip_on_cran()
  
  # have to pass in rfile
  expect_error(bigsearch(genus = "adfadf", email = "myrmecocystus@gmail.com", verbose=FALSE), 'argument "rfile" is missing')
  # have to pass in email
  expect_error(bigsearch(genus = "adfadf", rfile = "file1", verbose=FALSE), 'argument "email" is missing')
  # have to pass in at least one query term
  expect_error(bigsearch(rf = "file2", email = "stuff@stuff.com", verbose=FALSE))
})