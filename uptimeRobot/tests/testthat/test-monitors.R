test_that("uptimerobot.monitors", {
  
  skip_on_cran()
  
  api.key <- Sys.getenv("KEY", "")
  
  # Error in case the api.key is NULL or empty
  expect_error(uptimerobot.monitors("", "api.key cannot be empty or NULL"))
  
  # Error in case of invalid api.key
  expect_error(uptimerobot.monitors("blahblah"), "apiKey not mentioned or in a wrong format")
  
  # Get fields (typical)
  monitors.typical <- uptimerobot.monitors(api.key)
  
  # Output is a data frame
  expect_is(monitors.typical, "data.frame")
  
  # Output has the expected columns (typical)
  expect_identical(sort(colnames(monitors.typical)), sort(uptimerobot.fields("monitor")$typical))

  # Get fields (full)
  monitors.full <- uptimerobot.monitors(api.key, fields=uptimerobot.fields("monitor")$full)
  
  # Output is a data frame
  expect_is(monitors.full, "data.frame")
  
  # Output has the expected of columns (full)
  expect_identical(sort(colnames(monitors.full)), sort(uptimerobot.fields("monitor")$full))

  # Get fields (compact)
  monitors.compact <- uptimerobot.monitors(api.key, fields=uptimerobot.fields("monitor")$compact)
  
  # Output is a data frame
  expect_is(monitors.compact, "data.frame")
  
  # Output has the expected number of columns (typical)
  expect_identical(sort(colnames(monitors.compact)), sort(uptimerobot.fields("monitor")$compact))
  
  # Clean the environment
  rm(list = ls())
  
})
