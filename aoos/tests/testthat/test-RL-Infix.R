context("Infix")

test_that("Infix operators for S3", {
  
  Test <- function(.x) {
    ".+" <- function(e2) Test(getX() + e2)
    ".==" <- function(e2) getX() == e2
    ".>=" <- function(e2) getX() >= e2
    ".-" <- function(e2) {
      if (missing(e2)) -.x else .x - e2
    }
    getX <- function() .x
    retList(c("Test", "Infix"))
  }
  
  expect_equal(-Test(3), -3)
  expect_equal(Test(3) - 2, 1)
  expect_true(Test(2) + 2 == 4)
  expect_true(Test(2) == 2)
  expect_false(Test(3) == 2)
  expect_true(Test(2) >= 1)
  expect_true(Test(2) >= 2)
  expect_false(Test(2) >= 3)
  
})

test_that("Infix selects next method.", {
  Test <- function(.x) {
    ".+" <- function(e2) Test(getX() + e2)
    ".==" <- function(e2) getX() == e2
    getX <- function() .x
    retList(c("Test", "Infix"))
  }
  
  expect_true(Test(2) + 2 == 4)
  # This is not really testing what I want. It should test that S3 dispatch is 
  # continued after checking for infix functions encapsulated in class. Error
  # message depends on language settings...
  expect_error(Test(2) - 2 == 4)
  
})

