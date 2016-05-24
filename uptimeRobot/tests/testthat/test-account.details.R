test_that("uptimerobot.account.details", {
  
  skip_on_cran()
  
  api.key <- Sys.getenv("KEY", "")
  
  account.details.list <- uptimerobot.account.details(api.key)
  
  # Output is a list
  expect_is(account.details.list, "list")
  
  # Output list has 5 elements
  expect_equal(length(account.details.list), 5)
  
  # Elements are named as expected
  expect_identical(names(account.details.list), c("monitorLimit", "monitorInterval", "upMonitors", "downMonitors", "pausedMonitors"))
  
  account.details.vector <- uptimerobot.account.details(api.key, unlist = TRUE)
  
  # Output is a vector
  expect_is(account.details.vector, "integer")
    
  # Output vector has 5 rows
  expect_equal(length(account.details.vector), 5)
  
  # Error in case of invalid api.key
  expect_error(uptimerobot.account.details("blahblah"), "apiKey not mentioned or in a wrong format")
  
  # Error in case the api.key is NULL or empty
  expect_error(uptimerobot.account.details("", "api.key cannot be empty or NULL"))
  
  # Clean the environment
  rm(list = ls())
  
})
