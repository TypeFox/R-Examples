context('run_var')

test_that('run_var calls and returns the VAR function correctly', {
  called_count_run_var <<- 0
  called_count_restrictions <<- 0
  expected_result <<- list(datamat = matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, c('a', 'b', 'a1', 'b1', 'const'))),
                           y = matrix(NA, nrow = 1, ncol = 2, dimnames = list(NULL, c('a', 'b'))))
  with_mock(
    `vars::VAR` = function(...) {
      called_count_run_var <<- called_count_run_var + 1
      expect_equal(list(...), list(y = 1, p = 3, exogen = 2))
      expected_result
    },
    `autovarCore:::restrictions_for_lag` = function(...) {
      called_count_restrictions <<- called_count_restrictions + 1
      expect_equal(list(...), list(3, 2, 3))
      NULL
    },
    expect_equal(autovarCore:::run_var(1, 2, 3), expected_result)
  )
  expect_equal(called_count_run_var, 1)
  expect_equal(called_count_restrictions, 1)
  rm(list = c('called_count_run_var',
              'called_count_restrictions',
              'expected_result'), pos = '.GlobalEnv')
})

test_that('run_var calls the restrict function correctly for lag 0 models', {
  called_count_run_var <<- 0
  called_count_restrictions <<- 0
  called_count_restrict <<- 0
  expected_result <<- list(datamat = matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, c('a', 'b', 'a1', 'b1', 'const'))),
                           y = matrix(NA, nrow = 1, ncol = 2, dimnames = list(NULL, c('a', 'b'))),
                           varresult = list(list(terms = "a"), list(terms = "b")))
  expected_result2 <<- expected_result
  attr(expected_result2$varresult[[1]]$terms, "intercept") <- 1
  attr(expected_result2$varresult[[2]]$terms, "intercept") <- 1
  returned_resmat <<- matrix(NA, ncol = 3, nrow = 2)
  returned_resmat[1, ] <- c(0, 0, 1)
  returned_resmat[2, ] <- c(0, 0, 1)
  with_mock(
    `vars::VAR` = function(...) {
      called_count_run_var <<- called_count_run_var + 1
      expect_equal(list(...), list(y = 3, p = 1, exogen = 2))
      expected_result
    },
    `autovarCore:::restrictions_for_lag` = function(...) {
      called_count_restrictions <<- called_count_restrictions + 1
      expect_equal(list(...), list(0, 2, 3))
      returned_resmat
    },
    `vars::restrict` = function(...) {
      called_count_restrict <<- called_count_restrict + 1
      expect_equal(list(...), list(expected_result, method = "manual", resmat = returned_resmat))
      expected_result
    },
    expect_equal(autovarCore:::run_var(3, 2, 0), expected_result2)
  )
  expect_equal(called_count_run_var, 1)
  expect_equal(called_count_restrictions, 1)
  expect_equal(called_count_restrict, 1)
  rm(list = c('called_count_run_var',
              'called_count_restrictions',
              'expected_result',
              'expected_result2',
              'called_count_restrict',
              'returned_resmat'), pos = '.GlobalEnv')
})

test_that('run_var calls the restrict function correctly for lag 2 models', {
  called_count_run_var <<- 0
  called_count_restrict <<- 0
  expected_result <<- list(datamat = matrix(NA, nrow = 1, ncol = 12, dimnames = list(NULL, c('a', 'b', 'c', 'a1', 'b1', 'c1', 'a2', 'b2', 'c2', 'const', 'dailymeas_1', 'dailymeas_2'))),
                           y = matrix(NA, nrow = 1, ncol = 3, dimnames = list(NULL, c('a', 'b', 'c'))),
                           varresult = list(list(terms = "a"), list(terms = "b"), list(terms = "c")))
  expected_result2 <<- expected_result
  attr(expected_result2$varresult[[1]]$terms, "intercept") <- 1
  attr(expected_result2$varresult[[2]]$terms, "intercept") <- 1
  attr(expected_result2$varresult[[3]]$terms, "intercept") <- 1
  returned_resmat <<- matrix(NA, ncol = 9, nrow = 3)
  returned_resmat[1, ] <<- c(1, 1, 1, 1, 0, 0, 1, 1, 1)
  returned_resmat[2, ] <<- c(1, 1, 1, 0, 1, 0, 1, 1, 1)
  returned_resmat[3, ] <<- c(1, 1, 1, 0, 0, 1, 1, 1, 1)
  with_mock(
    `vars::VAR` = function(...) {
      called_count_run_var <<- called_count_run_var + 1
      expect_equal(list(...), list(y = 3, p = 2, exogen = 1))
      expected_result
    },
    `vars::restrict` = function(...) {
      called_count_restrict <<- called_count_restrict + 1
      expect_equal(list(...), list(expected_result, method = "manual", resmat = returned_resmat))
      expected_result
    },
    expect_equal(autovarCore:::run_var(3, 1, 2), expected_result2)
  )
  expect_equal(called_count_run_var, 1)
  expect_equal(called_count_restrict, 1)
  rm(list = c('called_count_run_var',
              'expected_result',
              'expected_result2',
              'called_count_restrict',
              'returned_resmat'), pos = '.GlobalEnv')
})

test_that('restrictions_for_lag returns no restrictions for lag 1 or 3+', {
  expect_null(autovarCore:::restrictions_for_lag(1, 3, 4))
  expect_null(autovarCore:::restrictions_for_lag(3, 3, 10))
  expect_null(autovarCore:::restrictions_for_lag(7, 2, 21))
})

test_that('restrictions_for_lag returns the correct restrictions for lag 0 models', {
  expected_result <- matrix(NA, ncol = 8, nrow = 2)
  expected_result[1, ] <- c(0, 0, 1, 1, 1, 1, 1, 1)
  expected_result[2, ] <- c(0, 0, 1, 1, 1, 1, 1, 1)
  expect_equal(autovarCore:::restrictions_for_lag(0, 2, 8), expected_result)
  expected_result <- matrix(NA, ncol = 5, nrow = 4)
  expected_result[1, ] <- c(0, 0, 0, 0, 1)
  expected_result[2, ] <- c(0, 0, 0, 0, 1)
  expected_result[3, ] <- c(0, 0, 0, 0, 1)
  expected_result[4, ] <- c(0, 0, 0, 0, 1)
  expect_equal(autovarCore:::restrictions_for_lag(0, 4, 5), expected_result)
})

test_that('restrictions_for_lag returns the correct restrictions for lag 2 models', {
  expected_result <- matrix(NA, ncol = 5, nrow = 2)
  expected_result[1, ] <- c(1, 1, 1, 0, 1)
  expected_result[2, ] <- c(1, 1, 0, 1, 1)
  expect_equal(autovarCore:::restrictions_for_lag(2, 2, 5), expected_result)
  expected_result <- matrix(NA, ncol = 17, nrow = 5)
  expected_result[1, ] <- c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1)
  expected_result[2, ] <- c(1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1)
  expected_result[3, ] <- c(1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1)
  expected_result[4, ] <- c(1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1)
  expected_result[5, ] <- c(1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1)
  expect_equal(autovarCore:::restrictions_for_lag(2, 5, 17), expected_result)
  expected_result <- matrix(NA, ncol = 9, nrow = 3)
  expected_result[1, ] <- c(1, 1, 1, 1, 0, 0, 1, 1, 1)
  expected_result[2, ] <- c(1, 1, 1, 0, 1, 0, 1, 1, 1)
  expected_result[3, ] <- c(1, 1, 1, 0, 0, 1, 1, 1, 1)
  expect_equal(autovarCore:::restrictions_for_lag(2, 3, 9), expected_result)
})
