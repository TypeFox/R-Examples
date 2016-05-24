test_that("uptimerobot.contacts", {
  
  skip_on_cran()
  
  api.key <- Sys.getenv("KEY", "")
  
  # Error in case the api.key is NULL or empty
  expect_error(uptimerobot.contacts("", "api.key cannot be empty or NULL"))
  
  # Error in case of invalid api.key
  expect_error(uptimerobot.contacts("blahblah"), "apiKey not mentioned or in a wrong format")
  
  # Get fields (typical)
  contacts.typical <- uptimerobot.contacts(api.key)
  
  # Output is a data frame
  expect_is(contacts.typical, "data.frame")
  
  # Output has the expected columns (typical)
  expect_identical(sort(colnames(contacts.typical)), sort(uptimerobot.fields("contact")$typical))

  # Get fields (full)
  contacts.full <- uptimerobot.contacts(api.key, fields=uptimerobot.fields("contact")$full)
  
  # Output is a data frame
  expect_is(contacts.full, "data.frame")
  
  # Output has the expected columns (full)
  expect_identical(sort(colnames(contacts.full)), sort(uptimerobot.fields("contact")$full))

  # Get fields (compact)
  contacts.compact <- uptimerobot.contacts(api.key, fields=uptimerobot.fields("contact")$compact)
  
  # Output is a data frame
  expect_is(contacts.compact, "data.frame")
  
  # Output has the expected columns (typical)
  expect_identical(sort(colnames(contacts.compact)), sort(uptimerobot.fields("contact")$compact))
  
  # Clean the environment
  rm(list = ls())
  
})
