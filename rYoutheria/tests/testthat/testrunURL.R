context("Test runURL")

# Most things are tested elsewhere

test_that("Testing errors", {
  
  expect_error(rYoutheria:::runURL(URL = NULL, type = 'm'),
               'URL is NULL')
  
  expect_error(rYoutheria:::runURL(URL = '', type = NULL),
               'type is NULL')  
  
})
