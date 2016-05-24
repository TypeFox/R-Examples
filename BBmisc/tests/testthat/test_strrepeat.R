context("strrepeat")

test_that("strrepeat", {
  expect_identical(strrepeat("x", 3), "xxx")
  expect_identical(strrepeat("x", 3, "y"), "xyxyx")
  expect_identical(strrepeat(c("x", "y"), 2), "xyxy")
})

