
library(testthat)
library(imputeMulti)

context("int- supDist works")

test_that("supDist errors and results", {
  set.seed(315)
  x1 <- rnorm(10)
  x2 <- rnorm(100)
  y  <- rnorm(100)
  
  expect_error(supDist(x1, y)) 
  expect_error(supDist(y, x1))
  
  expect_equal(supDist(x2,y), max(abs(x2-y)))
})


#--------------------------------------
context("int- marg_compare works")

test_that("errors work; return type is correct", {
  x <- 1:5; dim(x) <- c(1,5)
  x <- rbind(x,x)
  x2 <- 1:6; dim(x2) <- c(1,6)
  x2 <- rbind(x2,x2)
  x3 <- x2[,1:5]
  
  expect_error(marg_complete_compare(x, x2, FALSE))
  expect_error(marg_complete_compare(x, x2, TRUE))
  
  expect_true(is.list(marg_comp_compare(x, x3, FALSE)))
  expect_true(is.list(marg_comp_compare(x, x3, TRUE)))
  expect_equal(length(marg_comp_compare(x, x3, FALSE)), 2)
  expect_equal(length(marg_comp_compare(x, x3, TRUE)), 2)
})

test_that("marg_compare works correctly", {
  # set up
  set.seed(125)
  x1 <- factor(sample(1:5, size=10, replace= TRUE))
  x2 <- factor(sample(6:10, size=10, replace= TRUE))
  x3 <- factor(sample(11:15, size=10, replace= TRUE))
  x4 <- factor(sample(16:20, size=10, replace= TRUE))
  x5 <- factor(sample(21:26, size=10, replace= TRUE))
  
  dat2 <- dat <- c(x1, x2, x3, x4, x5)
  # insert missing values
  mis.ind <- sample(1:length(dat), size= 10, replace= FALSE)
  dat[mis.ind] <- NA
  rm(x1,x2,x3,x4,x5, mis.ind)
  dim(dat)<- dim(dat2) <- c(10, 5)
  
  expect_equal(unlist(marg_comp_compare(dat, dat, TRUE)), 1:10)
  expect_equal(unlist(marg_comp_compare(dat, dat, FALSE)), 1:10)
  expect_equal(unlist(marg_comp_compare(dat2, dat2, TRUE)), 1:10)
  expect_equal(unlist(marg_comp_compare(dat2, dat2, FALSE)), 1:10)
  expect_equal(unlist(marg_comp_compare(dat, dat2, TRUE)), 1:10)
  expect_equal(unlist(marg_comp_compare(dat, dat2, FALSE)), 1:10)
  
  set.seed(125)
  dat2 <- data.frame(x1= factor(sample(1:5, size=10, replace= TRUE)),
                     x2= factor(sample(6:10, size=10, replace= TRUE)),
                     x3= factor(sample(11:15, size=10, replace= TRUE)),
                     x4= factor(sample(16:20, size=10, replace= TRUE)),
                     x5= factor(sample(21:26, size=10, replace= TRUE)))
  
  # insert missing values
  dat <- dat2
  dat[c(3,7),1] <- NA
  dat[c(1,5),2] <- NA
  dat[c(2,6),3] <- NA
  dat[c(7,9),4] <- NA
  dat[c(4,10),5] <- NA
  
  expect_equal(unlist(marg_complete_compare(dat, dat, TRUE)), c(1:10,9,10))
  expect_equal(unlist(marg_complete_compare(dat, dat, FALSE)), c(1:10,9,10))
  expect_equal(unlist(marg_complete_compare(dat2, dat2, TRUE)), 1:10)
  expect_equal(unlist(marg_complete_compare(dat2, dat2, FALSE)), 1:10)
  expect_equal(unlist(marg_complete_compare(dat, dat2, TRUE)), 1:10)
  expect_equal(unlist(marg_complete_compare(dat, dat2, FALSE)), 1:10)
  
})


