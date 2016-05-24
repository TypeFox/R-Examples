context('select_valid_masks')

test_that('select_valid_masks returns the correct result', {
  all_outlier_masks <- 0:7
  invalid_mask <- 3 # columns 1 and 2 are invalid
  expected_result <- c(0, 4)
  expect_equal(autovarCore:::select_valid_masks(all_outlier_masks,
                                                invalid_mask),
               expected_result)
})

test_that('select_valid_masks calls bitwAnd correctly', {
  called_count <<- 0
  input_masks <<- 0:31
  with_mock(
    `base::bitwAnd` = function(...) {
      called_count <<- called_count + 1
      expect_equal(list(...), list(input_masks[called_count], 1))
      0
    },
    expect_equal(autovarCore:::select_valid_masks(input_masks, 1),
                 input_masks)
  )
  expect_equal(called_count, 32)
  rm(list = c('called_count', 'input_masks'), pos = '.GlobalEnv')
})

test_that('select_valid_masks works with an empty invalid_mask', {
  expect_equal(autovarCore:::select_valid_masks(0:7, 0), 0:7)
  expect_equal(autovarCore:::select_valid_masks(0:63, 0), 0:63)
})
