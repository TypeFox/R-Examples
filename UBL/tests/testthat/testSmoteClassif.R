context("SmoteClassif")

test_that("smoteClassif return a data frame", {
  data(iris)
  ir <- iris[c(1:70, 100:135), ]
  irnew <- SmoteClassif(Species~., ir)
  expect_true(is.data.frame(irnew))
})
