context("Test for lets.maplizer")

data(PAM)
data(IUCN)
trait <- IUCN$Description_Year

test_that("lets.maplizer works fine", {
  
  resu_test <- lets.maplizer(PAM, trait, PAM$S)
  expect_equal(class(resu_test), "matrix")
  expect_true(ncol(resu_test) == 3)
})

test_that("lets.maplizer works fine, other func", {
  
  resu_test <- lets.maplizer(PAM, trait, PAM$S, func = sd)
  expect_equal(class(resu_test), "matrix")
  expect_true(ncol(resu_test) == 3)
})

test_that("lets.maplizer works fine", {
  
  resu_test <- lets.maplizer(PAM, trait, PAM$S, ras = TRUE)
  expect_equal(class(resu_test), "list")
  expect_equal(class(resu_test[[1]]), "matrix")
  expect_true(inherits(resu_test[[2]], "RasterLayer"))
  
  expect_true(ncol(resu_test[[1]]) == 3)
})
