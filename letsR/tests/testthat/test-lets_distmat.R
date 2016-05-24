context("Test for lets.distmat")
data(PAM)
dimPAM <- summary(PAM)$Numberofcells

coords <- PAM[[1]][1:10, 1:2]


test_that("lets.distmat works fine, asdist = TRUE", {
  
  
  distPAM <- lets.distmat(PAM)   
  expect_true(class(distPAM) == "dist")
  expect_true(all(dim(as.matrix(distPAM)) == dimPAM))
    
})


test_that("lets.distmat works fine, asdist = FALSE", {
  
  
  distPAM <- lets.distmat(PAM, asdist = FALSE)   
  expect_true(class(distPAM) == "matrix")
  expect_true(all(dim(distPAM) == dimPAM))
  
})


test_that("lets.distmat works fine, miles = TRUE", {
  
  
  distPAM <- lets.distmat(PAM, miles = TRUE)   
  expect_true(class(distPAM) == "dist")
  expect_true(all(dim(as.matrix(distPAM)) == dimPAM))
  
})

test_that("lets.distmat works fine, xy as matrix ", {
  
  
  distPAM <- lets.distmat(coords, miles = TRUE)   
  expect_true(class(distPAM) == "dist")
  expect_true(all(dim(as.matrix(distPAM)) == nrow(coords)))
  
})
