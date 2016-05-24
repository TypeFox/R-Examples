test_that("uptimerobot.monitor.new/edit/reset/delete", {
  
  skip_on_cran()
  
  api.key <- Sys.getenv("KEY", "")

  # Error in case of invalid api.key
  expect_error(uptimerobot.monitor.new("blahblah", friendly.name="valid.name", url="https://gabrielebaldassarre.com", type="http"), "apiKey not mentioned or in a wrong format")
  expect_error(uptimerobot.monitor.edit("blahblah", id=1234, friendly.name="valid.name", url="https://gabrielebaldassarre.com"), "apiKey not mentioned or in a wrong format")
  expect_error(uptimerobot.monitor.delete("blahblah", id=1234), "apiKey not mentioned or in a wrong format")
    
  # Error in case the api.key is NULL or empty
  expect_error(uptimerobot.monitor.new("", friendly.name="valid.name", url="https://gabrielebaldassarre.com", type="http"), "api.key cannot be empty or NULL")
  expect_error(uptimerobot.monitor.edit("", friendly.name="valid.name", url="https://gabrielebaldassarre.com"), "api.key cannot be empty or NULL")
  expect_error(uptimerobot.monitor.delete("", 1234), "api.key cannot be empty or NULL")
  
  # monitorID doesn't exist
  expect_error(uptimerobot.monitor.delete(api.key, 1234), "monitorID doesn't exist")
  
  # If parameters are valid, output is numeric
  monitor.id <- uptimerobot.monitor.new(api.key, friendly.name="monitor.test", url="https://gabrielebaldassarre.com", type="http")
  expect_is(monitor.id, "numeric")
  expect_more_than(monitor.id, 0)
  
  # The just-created contact can be edited...
  expect_true(uptimerobot.monitor.edit(api.key, monitor.id, friendly.name="Open Analytics"))
  
  # ...its stats can be reset anytime...
  expect_true(uptimerobot.monitor.reset(api.key, monitor.id))
  
  # ...and it can be deleted
  expect_true(uptimerobot.monitor.delete(api.key, monitor.id))
  
  # Clean the environment
  rm(list = ls())
  
})
