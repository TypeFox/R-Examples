test_that("uptimerobot.monitor.contacts", {
  
  skip_on_cran()
  
  api.key <- Sys.getenv("KEY", "")
  
  # Error in case the api.key is NULL or empty
  expect_error(uptimerobot.monitor.contacts("", "api.key cannot be empty or NULL"))
  
  # Error in case of invalid api.key
  expect_error(uptimerobot.monitor.contacts("blahblah", "invalid.monitor"), "apiKey not mentioned or in a wrong format")
  
  contacts.df <- uptimerobot.monitor.contacts(api.key, monitors=uptimerobot.monitors(api.key)[1,1])
  
  # Output is a data frame
  expect_is(contacts.df, "data.frame")

  # Output has the expected columns
  expect_identical(colnames(contacts.df), c("monitor.id", "id", "type", "value", "threshold", "recurrence"))
  
  # Clean the environment
  rm(list = ls())
  
})
