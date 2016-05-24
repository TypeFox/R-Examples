test_that("test.is_null.na.returns_false", {
  expect_false(assertive.properties::is_null(NA))
})

test_that("test.is_null.nan.returns_false", {
  expect_false(assertive.properties::is_null(NaN))
})

test_that("test.is_null.null.returns_true", {
  expect_true(assertive.properties::is_null(NULL))
}) 

test_that("test.is_not_null.na.returns_true", {
  expect_true(assertive.properties::is_not_null(NA))
})

test_that("test.is_not_null.nan.returns_true", {
  expect_true(assertive.properties::is_not_null(NaN))
})

test_that("test.is_not_null.null.returns_false", {
  expect_false(assertive.properties::is_not_null(NULL))
}) 
