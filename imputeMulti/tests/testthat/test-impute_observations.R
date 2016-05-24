
library(testthat)
library(imputeMulti)

context("int- impute multinomial")

test_that("missing value imputation works", {
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
  
  
  # 03. Run 6 iterations, imputate missing values
  #------------------------------------
  iter6 <- multinomial_em(x_y= x_y, z_Os_y= z_Os_y, enum_comp= enum_comp, 
                          n_obs= nrow(dat), conj_prior= "none", 
                          verbose= FALSE, max_iter= 6) 
  
  dat_miss2 <- impute_multinomial_all(dat_miss, iter6@mle_x_y, p=ncol(dat))
  imputed_data <- rbind(dat_comp, dat_miss2)
  
  ### tests ###
  #------------------------------------
  expect_equal(sum(!complete.cases(imputed_data)), 0)
  expect_equal(sum(complete.cases(dat_miss2)), nrow(dat_miss2))
  for (i in 1:ncol(dat_miss2)) {
    expect_equal(levels(dat_miss2[,i]), levels(dat[,i]))
  }
  
})



