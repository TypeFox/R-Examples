context("Test for lets.subsetPAM")

data(PAM)

test_that("lets.subsetPAM works fine, remove.cells = TRUE", {
  
  
  resu_test <- lets.subsetPAM(PAM, PAM[[3]][1:20])
  
  expect_equal(class(resu_test), "PresenceAbsence")
  expect_equal(class(resu_test[[1]]), "matrix")
  expect_true(inherits(resu_test[[2]], "RasterLayer"))
  expect_equal(class(resu_test[[3]]), "character")
  
  response <- summary(resu_test)
  expect_true(response$Cellswithoutanypresence == 0)
  
})

test_that("lets.subsetPAM works fine, remove.cells = FALSE", {
  
  
  resu_test <- lets.subsetPAM(PAM, PAM[[3]][1:20],
                              remove.cells = FALSE)
  
  expect_equal(class(resu_test), "PresenceAbsence")
  expect_equal(class(resu_test[[1]]), "matrix")
  expect_true(inherits(resu_test[[2]], "RasterLayer"))
  expect_equal(class(resu_test[[3]]), "character")
  
  response <- summary(resu_test)
  expect_true(response$Cellswithoutanypresence > 0)
  
})

test_that("lets.subsetPAM error check", {
  
  expect_error(lets.subsetPAM(PAM, "Bruno"))
  expect_error(lets.subsetPAM(PAM[[1]], PAM[[3]][1:20]))
  expect_error(lets.subsetPAM(PAM, 1))
  expect_error(lets.subsetPAM(PAM[[1]], PAM[[3]][1:20]))
  expect_error(lets.subsetPAM(PAM, PAM[[3]][1:20], 1))  
  expect_warning(lets.subsetPAM(PAM, c("Bruno", PAM[[3]][1:20])))
  
})
