context("isScalarValue")

test_that("isScalarValue", {
  expect_true(isScalarValue(1))
  expect_true(isScalarValue(1L))
  expect_true(isScalarValue("a"))
  expect_true(isScalarValue(factor("a")))
  expect_true(isScalarValue(as.complex(1)))
  expect_true(isScalarValue(NA))
  
  expect_true(isScalarNumeric(1))
  expect_true(isScalarInteger(1L))
  expect_true(isScalarCharacter("a"))
  expect_true(isScalarFactor(factor("a")))
  expect_true(isScalarComplex(as.complex(1)))
  expect_true(isScalarLogical(NA))
  
  expect_false(isScalarComplex(1L))
  expect_false(isScalarInteger(1))
  expect_false(isScalarFactor("a"))
  expect_false(isScalarCharacter(factor("a")))
  expect_false(isScalarNumeric(as.complex(1)))
  expect_false(isScalarInteger(NA))
  
  expect_false(isScalarValue(NULL))
  expect_false(isScalarValue(iris))
  expect_false(isScalarValue(1:2))
  expect_false(isScalarValue(list(1)))
  
  expect_true(isScalarValue(NULL, null.ok=TRUE))
  expect_false(isScalarValue(NULL, na.ok=FALSE))
})
