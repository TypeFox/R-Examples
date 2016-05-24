context("padding_value")

test_that("padding_value can deal with all the possible inputs", {
  data("session_dataset")
  session_dataset$timestamp <- to_seconds(x = session_dataset$timestamp, format = "%Y%m%d%H%M%S")
  events_by_user <- split(session_dataset$timestamp, session_dataset$UUID)
  sessions <- reconstruct_sessions(events_by_user)
  
  expect_that(padding_value(sessions, "geometric mean"), equals(41.07547))
  expect_that(padding_value(sessions, "other", function(x,y){return(sum(x)+y)}, y = 19), equals(9656260))
  expect_that(padding_value(sessions, "arithmetic mean"),equals(226.115))
  expect_that(padding_value(sessions, "median"), equals(36))
})

test_that("padding_value can detect invalid inputs", {
  expect_that(padding_value("a turnip"), throws_error("VECTOR_ELT"))
  expect_that(padding_value(list(c(1,2,3)), "also a turnip"), throws_error("is missing, with no default"))
})