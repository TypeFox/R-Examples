context("sessionisation")

test_that("sessionisation handles zero-event entries", {
  
  #Example set.
  events <- list()
  
  #Test resulting value
  expect_that(length(reconstruct_sessions(events)), equals(0))
  
})

test_that("sessionisation handles single-event entries", {
  
  #Example set.
  events <- list(1)
  
  #Test resulting value
  expect_that(length(reconstruct_sessions(events)), equals(1))
  
})

test_that("sessionisation handles multi-set entries", {
  
  #Example set.
  events <- list(1,2)
  
  #Test resulting value
  expect_that(length(reconstruct_sessions(events)), equals(2))
  
})

test_that("sessionisation handles thresholds", {
  
  #Example set and run
  events <- list(c(1,2,3,12,5),2)
  results <- reconstruct_sessions(events,6)
  
  #Test resulting value
  expect_that(length(results), equals(3))
  expect_that(length(results[[1]]), equals(4))
  expect_that(length(results[[2]]), equals(1))
  expect_that(length(results[[3]]), equals(1))
  
  
})