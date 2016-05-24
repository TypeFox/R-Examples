test_that("uptimerobot.contact.new/delete", {
  
  skip_on_cran()
  
  api.key <- Sys.getenv("KEY", "")

  # Error in case of invalid api.key
  expect_error(uptimerobot.contact.new("blahblah", "email", value="jdoe@localdomain.it", "John Doe"), "apiKey not mentioned or in a wrong format")
  expect_error(uptimerobot.contact.delete("blahblah", 1234), "apiKey not mentioned or in a wrong format")
  
  # Error in case the api.key is NULL or empty
  expect_error(uptimerobot.contact.new(""), "api.key cannot be empty or NULL")
  expect_error(uptimerobot.contact.delete("", 1234), "api.key cannot be empty or NULL")
  
  # ContactID doesn't exist
  expect_error(uptimerobot.contact.delete(api.key, 1234), "alertContactID doesn't exist")
  
  # Error in case the type is empty or not recognized
  expect_error(uptimerobot.contact.new(api.key, ""), "contact type missing or not recognized.")
  expect_error(uptimerobot.contact.new(api.key, "blahblahblah"), "contact type missing or not recognized.")
  
  # Error in case the email is not valid
  expect_error(uptimerobot.contact.new(api.key, "email", value="invalid@localdomain", "John Doe"), "alertContactValue should be a valid e-mail for this alertContactType")
  
  # If parameters are valid, output is numeric
  contact.id <- uptimerobot.contact.new(api.key, "email", "jdoe@localdomain.it", "John Doe")
  expect_is(contact.id, "numeric")
  expect_more_than(contact.id, 0)
  
  # The just-created contact can be deleted
  expect_true(uptimerobot.contact.delete(api.key, contact.id))
  
  # Clean the environment
  rm(list = ls())
  
})
