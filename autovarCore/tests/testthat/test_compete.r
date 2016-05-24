context('compete')

test_that('compete calls its subfunctions correctly', {
  called_count <<- 0
  with_mock(
    `autovarCore:::challenger_wins` = function(...) {
      called_count <<- called_count + 1
      expect_equal(list(...), list(1, 2, 3))
      TRUE
    },
    expect_equal(autovarCore:::compete(1, 2, 3), 2)
  )
  expect_equal(called_count, 1)
  called_count <<- 0
  with_mock(
    `autovarCore:::challenger_wins` = function(...) {
      called_count <<- called_count + 1
      expect_equal(list(...), list(1, 2, 3))
      FALSE
    },
    expect_equal(autovarCore:::compete(1, 2, 3), 1)
  )
  expect_equal(called_count, 1)
  rm(list = 'called_count', pos = '.GlobalEnv')
})

test_that('challenger_wins calls its subfunctions correctly', {
  best <- list(bucket = 0.05, model_score = 100,
               nr_dummy_variables = 1)
  challenger <- list(bucket = 0.05, model_score = 99,
                     nr_dummy_variables = 2)
  expect_equal(autovarCore:::challenger_wins(best, challenger, TRUE), FALSE)
  best <- list(bucket = 0.05, model_score = 100,
               nr_dummy_variables = 1)
  challenger <- list(bucket = 0.05, model_score = 99,
                     nr_dummy_variables = 2)
  expect_equal(autovarCore:::challenger_wins(best, challenger, FALSE), TRUE)
})

test_that('challenger_wins returns the model with the highest bucket', {
  best <- list(bucket = 0.05)
  challenger <- list(bucket = 0.01)
  expect_equal(autovarCore:::challenger_wins(best, challenger, TRUE), FALSE)
})

test_that('challenger_wins returns the model with the least outlier columns if buckets are equal', {
  best <- list(bucket = 0.05, model_score = 99,
               nr_dummy_variables = 2)
  challenger <- list(bucket = 0.05, model_score = 100,
                     nr_dummy_variables = 1)
  expect_equal(autovarCore:::challenger_wins(best, challenger, TRUE), TRUE)
})

test_that('challenger_wins otherwise returns the model with lowest model_score', {
  best <- list(bucket = 0.05, model_score = 100,
               nr_dummy_variables = 1)
  challenger <- list(bucket = 0.05, model_score = 99,
               nr_dummy_variables = 1)
  expect_equal(autovarCore:::challenger_wins(best, challenger, TRUE), TRUE)
  expect_equal(autovarCore:::challenger_wins(best, challenger, FALSE), TRUE)
})

