# tests #########################
########## 2/3/4PL Model ########
################################################################
context("Jackknife")

# ------------------------- testing 1>>>


set.seed(1622)
# intercepts
diffpar <- seq(-3,3,length=12)
# slope parameters
sl     <- round(runif(12,0.5,1.5),2)
la     <- round(runif(12,0,0.25),2)
ua     <- round(runif(12,0.8,1),2)

# response matrix
awm <- matrix(sample(0:1,10*12,replace=TRUE),ncol=12)

awm <- rbind(awm,c(rep(1,11),0))
awm <- rbind(awm,rep(1,12))
awm <- rbind(awm,rep(1,12))



## 1PL model ##### 

# MLE estimation
res1plmle <- PP_4pl(respm = awm,thres = diffpar, slopes = rep(1,length(diffpar)),type = "mle")
# WLE estimation
res1plwle <- PP_4pl(respm = awm,thres = diffpar, slopes = rep(1,length(diffpar)),type = "wle")
# MAP estimation
res1plmap <- PP_4pl(respm = awm,thres = diffpar, slopes = rep(1,length(diffpar)),type = "map")



res_jk1 <- JKpp(res1plmle)
res_jk2 <- JKpp(res1plwle)
res_jk3 <- JKpp(res1plmap)


test_that("result = matrix with same dimension",{
  expect_that(dim(res1plmle$resPP$resPP),is_identical_to(dim(res_jk1$resjk)))
  expect_that(dim(res1plwle$resPP$resPP),is_identical_to(dim(res_jk2$resjk)))
  expect_that(dim(res1plmap$resPP$resPP),is_identical_to(dim(res_jk3$resjk)))
})





