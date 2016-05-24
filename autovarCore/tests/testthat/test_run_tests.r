context('run_tests')

test_that('run_tests returns the correct result', {
  called_count_run_test <<- 0
  called_count_portmanteau <<- 0
  test_names <<- c('a', 'b')
  with_mock(
    `autovarCore:::assess_portmanteau` = function(...) {
      called_count_portmanteau <<- called_count_portmanteau + 1
      expect_equal(list(...), list(1))
      called_count_portmanteau
    },
    `autovarCore:::run_test` = function(...) {
      called_count_run_test <<- called_count_run_test + 1
      expect_equal(list(...), list(test_names[[called_count_run_test]]))
      autovarCore:::assess_portmanteau
    },
    expect_equal(autovarCore:::run_tests(1, test_names), c(1, 2))
  )
  expect_equal(called_count_portmanteau, 2)
  expect_equal(called_count_run_test, 2)
  rm(list = c('called_count_portmanteau',
              'called_count_run_test',
              'test_names'), pos = '.GlobalEnv')
})

test_that('run_tests calls run_test correctly', {
  called_count_portmanteau <<- 0
  called_count_portmanteau_squared <<- 0
  called_count_skewness <<- 0
  called_count_kurtosis <<- 0
  test_names <- c('portmanteau', 'portmanteau_squared',
                  'skewness', 'kurtosis')
  with_mock(
    `autovarCore:::assess_portmanteau` = function(...) {
      called_count_portmanteau <<- called_count_portmanteau + 1
      expect_equal(list(...), list(1))
      1
    },
    `autovarCore:::assess_portmanteau_squared` = function(...) {
      called_count_portmanteau_squared <<- called_count_portmanteau_squared + 1
      expect_equal(list(...), list(1))
      2
    },
    `autovarCore:::assess_skewness` = function(...) {
      called_count_skewness <<- called_count_skewness + 1
      expect_equal(list(...), list(1))
      3
    },
    `autovarCore:::assess_kurtosis` = function(...) {
      called_count_kurtosis <<- called_count_kurtosis + 1
      expect_equal(list(...), list(1))
      4
    },
    expect_equal(autovarCore:::run_tests(1, test_names), c(1, 2, 3, 4))
  )
  expect_equal(called_count_portmanteau, 1)
  expect_equal(called_count_portmanteau_squared, 1)
  expect_equal(called_count_skewness, 1)
  expect_equal(called_count_kurtosis, 1)
  rm(list = c('called_count_portmanteau',
              'called_count_portmanteau_squared',
              'called_count_skewness',
              'called_count_kurtosis'), pos = '.GlobalEnv')
})
