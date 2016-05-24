

library(testthat)
library(imputeMulti)

context("multinomial DA")


test_that("basic error and CP error checks", {
  expect_error(multinomial_em(conj_prior= "foo"))
  expect_error(multinomial_em(conj_prior= c("none", "data.dep")))
  expect_error(multinomial_em(conj_prior= "data.dep", alpha= NULL))
  expect_error(multinomial_em(conj_prior= "flat.prior", alpha= NULL))
  expect_error(multinomial_em(conj_prior= "flat.prior", alpha= c(1,2,3)))
  expect_error(multinomial_em(conj_prior= "flat.prior", alpha= rnorm(100)))
  
  ## data dependent: first set up test
  #------------------------------------
  set.seed(12315)
  x1 <- factor(sample(1:5, size=100, replace= TRUE))
  x2 <- factor(sample(6:10, size=100, replace= TRUE))
  x3 <- factor(sample(11:15, size=100, replace= TRUE))
  x4 <- factor(sample(16:20, size=100, replace= TRUE))
  x5 <- factor(sample(21:26, size=100, replace= TRUE))
  
  dat <- c(x1, x2, x3, x4, x5)
  mis.ind <- sample(1:length(dat), size= 75, replace= FALSE)
  dat[mis.ind] <- NA
  dim(dat)<- c(100, 5)
  rm(x1,x2,x3,x4,x5, mis.ind)
  
  dat <- data.frame(apply(dat, 2, function(x) as.factor(x)))
  enum <- expand.grid(sapply(dat, levels))
  
  ## now run test:
  expect_error(multinomial_data_aug(enum_comp= enum, conj_prior= "data.dep",
                              alpha= 1))
  expect_error(multinomial_data_aug(enum_comp= enum, conj_prior= "data.dep",
                              alpha= vector("numeric", length= nrow(enum))))
})


test_that("multinomial DA is converging", {
  ### set up testing inputs:
  # 01. data
  #------------------------------------
  set.seed(12315)
  x1 <- factor(sample(1:5, size=100, replace= TRUE))
  x2 <- factor(sample(6:10, size=100, replace= TRUE))
  x3 <- factor(sample(11:15, size=100, replace= TRUE))
  x4 <- factor(sample(16:20, size=100, replace= TRUE))
  x5 <- factor(sample(21:26, size=100, replace= TRUE))

  dat <- c(x1, x2, x3, x4, x5)
  mis.ind <- sample(1:length(dat), size= 75, replace= FALSE)
  dat[mis.ind] <- NA
  dim(dat)<- c(100, 5)
  rm(x1,x2,x3,x4,x5, mis.ind)

  dat <- data.frame(apply(dat, 2, function(x) as.factor(x)))
  # 02. sufficient statistics
  #------------------------------------
  enum <- expand.grid(sapply(dat, function(x) return(c(levels(x), NA))))
  enum_comp <- enum[complete.cases(enum),]
  enum_miss <- enum[!complete.cases(enum),]
  enum_miss <- enum_miss[apply(enum_miss, 1, function(x) !all(is.na(x))),] # not all missing
  rownames(enum_comp) <- 1:nrow(enum_comp) # y \in Y

  dat_comp <- dat[complete.cases(dat),]
  dat_miss <- dat[!complete.cases(dat),]
  # complete data sufficient statistics
  x_y     <- count_levels(dat_comp, enum_list= enum_comp, hasNA= "no")
  # missing data marginal sufficient statistics
  z_Os_y  <- count_levels(dat_miss, enum_list= enum_miss, hasNA= "count.miss")


  # 03. Run 20 iterations, make sure log-lik is increasing
  #------------------------------------
  iter1 <- multinomial_data_aug(x_y= x_y, z_Os_y= z_Os_y, enum_comp= enum_comp,
                          conj_prior= "none", verbose= FALSE, burnin= 1)
  iter20 <- multinomial_data_aug(x_y= x_y, z_Os_y= z_Os_y, enum_comp= enum_comp,
                          conj_prior= "none", verbose= FALSE, burnin= 20)
  iter10 <- multinomial_data_aug(x_y= x_y, z_Os_y= z_Os_y, enum_comp= enum_comp,
                          conj_prior= "none", verbose= FALSE, burnin= 10)

  # tests:
  expect_lt(iter1@mle_log_lik, iter20@mle_log_lik)
  expect_lt(iter1@mle_log_lik, iter10@mle_log_lik)
  expect_null(iter20@mle_x_y$alpha) # no prior
  expect_null(iter20@mle_x_y$counts)
  expect_null(iter20@mle_x_y$theta_y1)
  expect_equal(iter20@mle_cp, "none")
  expect_equal(iter20@mle_iter, 20)
  expect_equal(sum(iter20@mle_x_y$theta_y), 1)

})