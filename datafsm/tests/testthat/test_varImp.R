library(datafsm)
context("Variable importance function")

test_that("var_imp() returns correct type of object", {
  state_mat <-  matrix(c(1,1,2,2, 1, 1, 2, 2),
                       nrow = 2, ncol = 2, byrow = FALSE)
  action_vec <- c(1,2)
  cdata <- matrix(c(NA,1,0,1,1, NA,0,1,0,0), nrow = 5, ncol = 2)
  outcome <- c(1,0,1,1,1)
  period <- 1:5
  res <- var_imp(state_mat, action_vec, cdata, outcome, period,
                 measure = "accuracy")
  expect_is(res, "numeric")
})

test_that("var_imp() returns numerically correct output and errors", {
  state_mat <-  matrix(c(1,1,2,2, 1, 1, 2, 2),
                       nrow = 2, ncol = 2, byrow = FALSE)
  action_vec <- c(1,2)
  cdata <- matrix(c(NA,1,0,1,1, NA,0,1,0,0), nrow = 5, ncol = 2)
  outcome <- c(1,0,1,1,1)
  period <- 1:5
  res <- var_imp(state_mat, action_vec, cdata, outcome, period,
                 measure = "accuracy")
  expect_equal(res[1], 100)
  
  state_mat <-  matrix(c(1,1,2,2, 1, 1, 2, 2),
                       nrow = 2, ncol = 2, byrow = FALSE)
  action_vec <- c(1,2)
  # cdata should cause error bc there is a row with more than one 1 in it
  cdata <- matrix(c(NA,1,0,1,1, NA,1,1,0,0), nrow = 5, ncol = 2)
  outcome <- c(1,0,1,1,1)
  period <- 1:5
  expect_error(var_imp(state_mat, action_vec, cdata, outcome, period,
                       measure = "accuracy"))
  
})

