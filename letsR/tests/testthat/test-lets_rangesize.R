context("Test for lets.rangesize")

data(PAM)
data(Phyllomedusa)
projection(Phyllomedusa) <- projection(PAM[[2]])

test_that("lets.rangesize works fine, geographic", {
  
  
  resu_test <- lets.rangesize(x = Phyllomedusa,
                              coordinates = "geographic")
  
  expect_equal(class(resu_test), "matrix")
  expect_true(ncol(resu_test) == 1)
  expect_true(!any(is.na(resu_test)))
})


test_that("lets.rangesize works fine, planar", {
  
  
  resu_test <- lets.rangesize(x = Phyllomedusa,
                              coordinates = "planar")
  
  expect_equal(class(resu_test), "matrix")
  expect_true(ncol(resu_test) == 1)
  expect_true(!any(is.na(resu_test)))
})

test_that("lets.rangesize works fine, squaremeter", {
  
  
  resu_test <- lets.rangesize(x = PAM,
                              units = "squaremeter")
  
  expect_equal(class(resu_test), "matrix")
  expect_true(ncol(resu_test) == 1)
  expect_true(!any(is.na(resu_test)))
})


test_that("lets.rangesize works fine, squaremeter", {
  
  
  resu_test <- lets.rangesize(x = PAM,
                              units = "cell")
  
  expect_equal(class(resu_test), "matrix")
  expect_true(ncol(resu_test) == 1)
  expect_true(!any(is.na(resu_test)))
})

test_that("lets.rangesize works fine, squaremeter global", {
  
  PAM2 <- lets.presab(Phyllomedusa, res = 5)
  resu_test <- lets.rangesize(x = PAM2,
                              units = "squaremeter")
  
  expect_equal(class(resu_test), "matrix")
  expect_true(ncol(resu_test) == 1)
  expect_true(!any(is.na(resu_test)))
})




test_that("lets.rangesize error", {
  
  
  expect_error(lets.rangesize(x = PAM, units = "cells"))
  
})

test_that("lets.rangesize error2", {
  
  
  expect_error(lets.rangesize(x = Phyllomedusa, coordinates = "cells"))
  
})
