context("session_events")

test_that("session_events handles zero-event entries", {
  
  #Example set.
  events <- list()
  
  #Test resulting value
  expect_that(length(session_events(events)), equals(0))
  
})

test_that("session_events handles single-event entries", {
  
  #Example set.
  events <- list(1)
  results <- session_events(events)
  
  #Test resulting value
  expect_that(length(results), equals(1))
  expect_that(results, equals(1))
  
})

test_that("session_events handles multi-set entries", {
  
  #Example set.
  events <- list(1,2)
  results <- session_events(events)
  
  #Test resulting value
  expect_that(length(results), equals(2))
  expect_that(results[1], equals(1))
  expect_that(results[2], equals(1))
  
})

test_that("session_events handles multi-event entries", {
  
  #Example set and run
  events <- list(c(1,2,3,12,5),2)
  results <- session_events(events)
  
  #Test resulting value
  expect_that(length(results), equals(2))
  expect_that(results[1], equals(5))
  expect_that(results[2], equals(1))
  
})