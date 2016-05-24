test_that("uptimerobot.monitor.logs", {
  
  skip_on_cran()
  
  api.key <- Sys.getenv("KEY", "")
  
  # Error in case the api.key is NULL or empty
  expect_error(uptimerobot.monitor.logs("", "api.key cannot be empty or NULL"))
  
  # Error in case of invalid api.key
  expect_error(uptimerobot.monitor.logs("blahblah", "invalid.monitor"), "apiKey not mentioned or in a wrong format")
  
  logs.df <- uptimerobot.monitor.logs(api.key, monitors=uptimerobot.monitors(api.key)[1,1])
  
  # Output is a data frame
  expect_is(logs.df, "data.frame")
  
  # Check if datetime is valid
  expect_equal(length(which(is.na(logs.df$datetime))), 0)
  
  # Output has the expected columns
  expect_identical(colnames(logs.df), c("monitor.id", "type", "datetime"))
  
  # Clean the environment
  rm(list = ls())
  
})
