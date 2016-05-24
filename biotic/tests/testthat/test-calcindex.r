library(biotic)
context("calcindex")

test_that("calcindex returns a dataframe", {
  # load built-in dataset
  data(almond)

  # calculate the BMWP index (default) for this dataset
  testindex<-calcindex(almond)
  # check that it outputs a dataframe object
  expect_is(testindex, "data.frame")

})

test_that("output rows equal input columns-1", {
  # load built-in dataset
  data(almond)

  # calculate the BMWP index (default) for this dataset
  testindex<-calcindex(almond)
  # check that it outputs a dataframe object
  expect_equal(nrow(testindex), ncol(almond)-1)

})

test_that("incorrect index specification causes an error", {
  # load built-in dataset
  data(almond)

  # calculate the BMWP index (default) for this dataset
  expect_error(calcindex(almond, index="FOO"))

})

test_that("abundance-weighted index causes warning with p-a data", {
  # create a dataset with just presence-absence data
  data(almond)
  presences<-almond[,-1]
  presences[presences>0]<-1
  p_a<-cbind.data.frame(almond[,1], presences)
  # calculate the LIFE index for this dataset (needs abundance data)
  expect_warning(calcindex(p_a, index="LIFE"))

})

test_that("PSI index for SeafieldDS1 = 62.50", {
  # load built-in dataset
  data(almond)

  # calculate the PSI index for this dataset
  testindex<-calcindex(almond, "PSI")
  # check that the value of PSI for SeafieldDS1 is 62.50
  expect_equal(testindex$PSI[5], 62.50)

})
