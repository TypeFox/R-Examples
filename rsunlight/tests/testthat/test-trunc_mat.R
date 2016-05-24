context("trunc_mat")

test_that("trunc_mat works with built in datasets", {
  skip_on_cran()

  iris_mod <- capture.output(trunc_mat(iris))

  # classes
  expect_is(iris, "data.frame")
  expect_is(iris_mod, "character")

  # values
  expect_equal(NROW(iris), 150)
  expect_equal(length(iris_mod), 12)
})

test_that("trunc_mat n parameter works", {
  skip_on_cran()

  iris_mod2 <- capture.output(trunc_mat(iris, n = 2))
  iris_mod100 <- capture.output(trunc_mat(iris, n = 100))

  expect_is(iris_mod2, "character")
  expect_is(iris_mod100, "character")
  expect_equal(length(iris_mod2), 4)
  expect_equal(length(iris_mod100), 102)
})
