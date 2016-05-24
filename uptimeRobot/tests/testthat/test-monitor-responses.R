test_that("uptimerobot.monitor.responses", {
  
  skip_on_cran()
  
  api.key <- Sys.getenv("KEY", "")
  
  # Error in case the api.key is NULL or empty
  expect_error(uptimerobot.monitor.responses("", "api.key cannot be empty or NULL"))
  
  # Error in case of invalid api.key
  expect_error(uptimerobot.monitor.responses("blahblah", "invalid.monitor"), "apiKey not mentioned or in a wrong format")
  
  responses.df <- uptimerobot.monitor.responses(api.key, monitors=uptimerobot.monitors(api.key)[1,1])
  
  # Output is a data frame
  expect_is(responses.df, "data.frame")
  
  # Check if datetime is valid
  expect_equal(length(which(is.na(responses.df$datetime))), 0)
  
  # Output has the expected columns
  expect_identical(colnames(responses.df), c("monitor.id", "datetime", "value"))
  
  # Clean the environment
  rm(list = ls())
  
})
