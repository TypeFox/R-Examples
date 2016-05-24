library(DNAprofiles)
context("KI distribution")
source("./ki_ibs_joint_dist_debug.R")

test_that(desc = "Compare the KI distribution with pencil and paper calculations",{
  
  set.seed(100)
  fr.true <- runif(10)
  fr.true <- fr.true/sum(fr.true)
  fr.ki <- runif(10)
  fr.ki <- fr.ki/sum(fr.ki)
    
  # homozygous case
  a <- 3
  b <- 3
  
  d1 <- Zdebug.ki.ibs.dist(a = a,b = b,hyp.1 = "FS",hyp.2 = "UN",hyp.true = "FS",f.ki = fr.ki,f.true = fr.true,
                     theta.ki = 0.03,theta.true = 0.03)

  d2 <- Zcond.ki.ibs.joint.dist.marker(a = a,b = b,k.hyp.1 = ibdprobs("FS"),k.hyp.2 = ibdprobs("UN"),k.hyp.true = ibdprobs("FS"),fr.ki = fr.ki,fr.true = fr.true,
                                 theta.ki = 0.03, theta.true = 0.03)
  
  
  # reorder and compare
  o <- sapply(d1[,1],function(y) which.min(abs(d2[,1]-y)))
  d2 <- d2[o,]  
  attr(d2,"dimnames") <- NULL
  expect_equal(d1,d2)
  
  # heterozygous case
  a <- 2
  b <- 9
  
  d1 <- Zdebug.ki.ibs.dist(a = a,b = b,hyp.1 = "FS",hyp.2 = "UN",hyp.true = "FS",f.ki = fr.ki,f.true = fr.true,
                                   theta.ki = 0.01,theta.true = 0.03)
  
  d2 <- Zcond.ki.ibs.joint.dist.marker(a = a,b = b,k.hyp.1 = ibdprobs("FS"),k.hyp.2 = ibdprobs("UN"),k.hyp.true = ibdprobs("FS"),fr.ki = fr.ki,fr.true = fr.true,
                                       theta.ki = 0.01, theta.true = 0.03)
    
  # reorder and compare
  o <- sapply(d1[,1],function(y) which.min(abs(d2[,1]-y)))
  d2 <- d2[o,]
  attr(d2,"dimnames") <- NULL
  
  expect_equal(d1,d2)  
})