context("collapsef")

test_that("collapsef", {
  expect_equal(collapsef("%s=%s", 1:2, 3:4), "1=3,2=4")
})
