context("Test getMeasurementTypes")

test_that("Testing error is given", {
  skip_on_cran()
  
  expect_error(getMeasurementTypes(measurementType = TRUE),
               regexp = 'argument must be numeric, integer, charater or NULL')
  expect_warning(getMeasurementTypes(measurementType = 'foo'),
               regexp = 'no matching measurement types found')
})

test_that("Testing data is returned", {
  skip_on_cran()
  
  expect_is(test <- getMeasurementTypes(), 'data.frame')
  expect_equal(ncol(test), 2)
  expect_true(all(c('Activity Cycle', 'Diet','Litter Size','Ranging Behaviour','Wing Morphology')
                  %in% test$Name))
  
  expect_is(test <- getMeasurementTypes(1), 'data.frame')
  expect_equal(ncol(test), 2)
  expect_equal(nrow(test), 1)
  expect_equal(test, data.frame(Id = 1, Name = 'Body Mass', stringsAsFactors = FALSE))
  
  expect_is(test <- getMeasurementTypes('Diet'), 'data.frame')
  expect_equal(ncol(test), 2)
  expect_equal(nrow(test), 1)
  expect_equal(test, data.frame(Id = 21, Name = 'Diet', stringsAsFactors = FALSE))
})