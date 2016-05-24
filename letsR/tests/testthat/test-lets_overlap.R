context("Test for lets.overlap")

data(PAM)


test_that("lets.overlap works fine, Chesser&Zink", {
  
  
  resu_test <- lets.overlap(PAM, method = "Chesser&Zink")  
  
  expect_equal(class(resu_test), "matrix")  
  expect_true(all(dim(resu_test) == length(PAM[[3]])))
  expect_true(!any(is.na(resu_test)))
  
})

test_that("lets.overlap works fine, xy = FALSE", {
  
  
  resu_test <- lets.overlap(PAM[[1]][, -c(1, 2)],
                            method = "Chesser&Zink",
                            xy = FALSE)  
  
  expect_equal(class(resu_test), "matrix")  
  expect_true(all(dim(resu_test) == length(PAM[[3]])))
  expect_true(!any(is.na(resu_test)))
})


test_that("lets.overlap works fine, xy = TRUE", {
  
  
  resu_test <- lets.overlap(PAM[[1]],
                            method = "Chesser&Zink",
                            xy = TRUE)  
  
  expect_equal(class(resu_test), "matrix")  
  expect_true(all(dim(resu_test) == length(PAM[[3]])))
  expect_true(!any(is.na(resu_test)))
})

test_that("lets.overlap works fine, Cells", {
  
  
  resu_test <- lets.overlap(PAM, method = "Cells")  
  
  expect_equal(class(resu_test), "matrix")  
  expect_true(all(dim(resu_test) == length(PAM[[3]])))
  expect_true(!any(is.na(resu_test)))
})

test_that("lets.overlap works fine, Proportional", {
  
  
  resu_test <- lets.overlap(PAM, method = "Proportional")  
  
  expect_equal(class(resu_test), "matrix")  
  expect_true(all(dim(resu_test) == length(PAM[[3]])))
  expect_true(!any(is.na(resu_test)))
})


test_that("lets.overlap error method", {
  
  expect_error(lets.overlap(PAM, method = "testerror"))
})

test_that("lets.overlap error xy", {
  
  expect_error(lets.overlap(PAM[[1]], method = "Cells"))
})