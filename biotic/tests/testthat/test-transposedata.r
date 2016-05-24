library(biotic)
context("transposedata")

test_that("transposedata returns a dataframe", {
  # load built-in dataset
  data(almond)

  # transpose this dataset
  t_almond<-transposedata(almond)
  # check that it outputs a dataframe object
  expect_is(t_almond, "data.frame")

})

test_that("transposed dataset rows equal input columns-1", {
  # load built-in dataset
  data(almond)

  # transpose this dataset
  t_almond<-transposedata(almond)
  # expect transposed rows = input cols -1
  expect_equal(nrow(t_almond), ncol(almond)-1)

})

