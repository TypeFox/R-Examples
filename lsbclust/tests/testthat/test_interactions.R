library("lsbclust")

## Check whether extended-precision long doubles are available
donttest <- .Machine$sizeof.longdouble == 0

## Check int.lsbclust output for consistency for different options of the "fixed" parameter

context("Interactions: fixed = 'columns'")

## Code to produce int.lsbclust object for fixed = "columns"
# data("dcars")
# set.seed(1)
# intcols <- int.lsbclust(data = dcars, margin = 3, nclust = 5, delta = c(1, 1, 1, 1), ndim = 2, 
#                            fixed = "columns", nstart = 2)
# save(intcols, file = "./tests/testthat/intcols.rda")

## Load saved object containing object 'intcols'
load("intcols.rda")

## Rerun code
data("dcars")
set.seed(1)
intcols.test <- int.lsbclust(data = dcars, margin = 3, nclust = 5, delta = c(1, 1, 1, 1), ndim = 2, 
                           fixed = "columns", nstart = 2, verbose = -1)

## Check different aspects of results
test_that("C matrices", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(intcols.test$C, intcols$C)
  })
test_that("D matrices", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(intcols.test$D, intcols$D)
  })
test_that("row fit", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(intcols.test$rfit, intcols$rfit)
  })
test_that("row fit per dimension", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(intcols.test$rfit.comp, intcols$rfit.comp)
  })
test_that("column fit", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(intcols.test$cfit, intcols$cfit)
  })
test_that("column fit per dimension", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(intcols.test$cfit.comp, intcols$cfit.comp)
  })
test_that("overall fit", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(intcols.test$ofit, intcols$ofit)
  })
test_that("cluster means", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(intcols.test$means, intcols$means)
  })
test_that("loss values", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(intcols.test$loss, intcols$loss)
  })
test_that("cluster membership", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_identical(intcols.test$cluster, intcols$cluster)
  })

## Remove loaded object
rm(intcols, intcols.test)

context("Interactions: fixed = 'rows'")

## Code to produce int.lsbclust object for fixed = "rows"
# data("dcars")
# set.seed(1)
# introws <- int.lsbclust(data = dcars, margin = 3, nclust = 8, delta = c(1, 0, 1, 1), ndim = 4, 
#                            fixed = "rows", nstart = 2)
# save(introws, file = "./tests/testthat/introws.rda")

## Load saved object containing object 'intcols'
load("introws.rda")

## Rerun code
data("dcars")
set.seed(1)
introws.test <- int.lsbclust(data = dcars, margin = 3, nclust = 8, delta = c(1, 0, 1, 1), 
                                    ndim = 4, fixed = "rows", nstart = 2, verbose = -1)

## Check different aspects of results
test_that("C matrices", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(introws.test$C, introws$C)
  })
test_that("D matrices", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(introws.test$D, introws$D)
  })
test_that("row fit", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(introws.test$rfit, introws$rfit)
  })
test_that("row fit per dimension", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(introws.test$rfit.comp, introws$rfit.comp)
  })
test_that("column fit", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(introws.test$cfit, introws$cfit)
  })
test_that("column fit per dimension", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(introws.test$cfit.comp, introws$cfit.comp)
  })
test_that("overall fit", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(introws.test$ofit, introws$ofit)
  })
test_that("cluster means", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(introws.test$means, introws$means)
  })
test_that("loss values", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(introws.test$loss, introws$loss)
  })
test_that("cluster membership", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_identical(introws.test$cluster, introws$cluster)
  })

## Remove loaded object
rm(introws, introws.test)

context("Interactions: fixed = 'none'")

## Code to produce int.lsbclust object for fixed = "columns"
# data("dcars")
# set.seed(12)
# intnone <- int.lsbclust(data = dcars, margin = 3, nclust = 6, delta = c(0, 1, 1, 1), ndim = 4, 
#                            fixed = "none", nstart = 2)
# save(intnone, file = "./tests/testthat/intnone.rda")

## Load saved object containing object 'intcols'
load("intnone.rda")

## Rerun code
data("dcars")
set.seed(12)
intnone.test <- int.lsbclust(data = dcars, margin = 3, nclust = 6, delta = c(0, 1, 1, 1), 
                                    ndim = 4, fixed = "none", nstart = 2, verbose = -1)

## Check different aspects of results
test_that("C matrices", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(intnone.test$C, intnone$C)
  })
test_that("D matrices", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(intnone.test$D, intnone$D)
  })
test_that("row fit", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(intnone.test$rfit, intnone$rfit)
  })
test_that("row fit per dimension", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(intnone.test$rfit.comp, intnone$rfit.comp)
  })
test_that("column fit", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(intnone.test$cfit, intnone$cfit)
  })
test_that("column fit per dimension", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(intnone.test$cfit.comp, intnone$cfit.comp)
  })
test_that("overall fit", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(intnone.test$ofit, intnone$ofit)
  })
test_that("cluster means", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(intnone.test$means, intnone$means)
  })
test_that("loss values", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_equal(intnone.test$loss, intnone$loss)
  })
test_that("cluster membership", {
  if (donttest) skip(".Machine$sizeof.longdouble is 0.")
  expect_identical(intnone.test$cluster, intnone$cluster)
  })

## Remove loaded objects
rm(intnone, intnone.test)