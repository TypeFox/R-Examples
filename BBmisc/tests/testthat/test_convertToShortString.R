context("convertToShortString")

test_that("convertToShortString", {
  expect_equal(convertToShortString(1L), "1")
  expect_equal(convertToShortString(1.0), "1")
  expect_equal(convertToShortString(1.23), "1.23")
  expect_equal(convertToShortString(numeric(0)), "numeric(0)")
  expect_equal(convertToShortString(factor(c())), "factor(0)")
  expect_equal(convertToShortString(iris), "<data.frame>")
  expect_equal(convertToShortString(list(a=1, 45)), "a=1, <unnamed>=45")
  expect_equal(convertToShortString(list(a=1, b=list(x=3))), "a=1, b=<list>")
  expect_equal(convertToShortString(list(a=1, b=iris)), "a=1, b=<data.frame>")
            
  expect_equal(convertToShortString(list()), "")
  expect_equal(convertToShortString(list(a=1)), "a=1")
  expect_equal(convertToShortString(list(a=1:2)), "a=1,2")
  expect_equal(convertToShortString(list(a=1:20)), "a=1,2,3,4,5,6,...")
  expect_equal(convertToShortString(list(a=1, 2, b=3)), "a=1, <unnamed>=2, b=3")
  expect_equal(convertToShortString(list(a=1, 2, b=data.frame())), "a=1, <unnamed>=2, b=<data.frame>")
  expect_equal(convertToShortString(list(a=identity, b=new.env())), "a=<function>, b=<environment>")

  expect_equal(convertToShortString(list(a=1, b=3.2)), "a=1, b=3.2")
  expect_equal(convertToShortString(list(a=1, b=3.223), num.format="%.2f"), "a=1.00, b=3.22")
  expect_equal(convertToShortString(list(a=1L, b=3.223), num.format="%.2f"), "a=1, b=3.22")
})


