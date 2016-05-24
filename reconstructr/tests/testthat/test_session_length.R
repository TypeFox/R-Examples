context("session_length")

test_that("session_length handles zero-event entries", {
  
  #Example set.
  events <- list()
  
  #Test resulting value
  expect_that(length(session_length(events)), equals(0))
  
})

test_that("session_length handles one-event entries", {
  
  #Example set.
  events <- list(1)
  result <- session_length(events)
  
  #Test resulting value
  expect_that(length(result), equals(1))
  expect_that(result[1], equals(-1))
  
  
})

test_that("session_length handles multi-event entries", {
  
  #Example set.
  events <- list(c(1,12))
  result <- session_length(events)
  
  #Test resulting value
  expect_that(length(result), equals(1))
  expect_that(result[1], equals(441))
  
  
})

test_that("session_length handles preserve_single_events correctly", {
  
  #Example set.
  events <- list(c(1))
  result <- session_length(events, preserve_single_events = TRUE)
  
  #Test resulting value
  expect_that(length(result), equals(1))
  expect_that(result[1], equals(430))
  
})

test_that("session_length handles strip_last correctly", {
  
  #Example set.
  events <- list(c(1))
  result <- session_length(events, strip_last = TRUE)
  
  #Test resulting value
  expect_that(length(result), equals(1))
  expect_that(result[1], equals(-1))
  
})

test_that("session_length handles variance in padding_values correctly", {
  
  #Example set.
  events <- list(c(1,12))
  result <- session_length(events, padding_value = 0)
  
  #Test resulting value
  expect_that(length(result), equals(1))
  expect_that(result[1], equals(11))
  
})

test_that("session_length handles preserve_single_events=F and strip_last=T correctly
          when dealing with two-event sessions.", {
  
  #Example set.
  events <- list(c(1,12))
  result <- session_length(events, strip_last = TRUE)
  
  #Test resulting value
  expect_that(length(result), equals(1))
  expect_that(result[1], equals(-1))
  
})