context("Test for lets.summarizer")

data(PAM)
data(temp)
pamvar <- lets.addvar(PAM, temp)


test_that("lets.summarizer works fine", {
  
  
  resu_test <- lets.summarizer(x = pamvar, pos = ncol(pamvar), xy = TRUE)

  expect_equal(class(resu_test), "data.frame")
  expect_true(nrow(resu_test) == length(PAM[[3]]))
  expect_true(ncol(resu_test) == 2)
})


test_that("lets.summarizer works fine, xy = FALSE", {
  
  
  resu_test <- lets.summarizer(x = pamvar[, -(1:2)], 
                               pos = ncol(pamvar[, -(1:2)]),
                               xy = FALSE)
  
  expect_equal(class(resu_test), "data.frame")
  expect_true(nrow(resu_test) == length(PAM[[3]]))
  expect_true(ncol(resu_test) == 2)
})


test_that("lets.summarizer works fine, other fun", {
  
  
  resu_test <- lets.summarizer(x = pamvar, pos = ncol(pamvar), 
                               xy = TRUE, fun = sd)
  
  expect_equal(class(resu_test), "data.frame")
  expect_true(nrow(resu_test) == length(PAM[[3]]))
  expect_true(ncol(resu_test) == 2)
})
