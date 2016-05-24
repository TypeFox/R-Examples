context("getClass1")

test_that("getClass1", {
  expect_equal(getClass1(iris), "data.frame")
  expect_equal(getClass1(1), "numeric")
  expect_equal(getClass1(NULL), "NULL")
  x = makeS3Obj(c("C1", "C2"), foo = 2)
  expect_equal(getClass1(x), "C1")
})


