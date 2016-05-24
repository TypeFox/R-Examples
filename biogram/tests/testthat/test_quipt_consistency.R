context("Exact test and distribution under H0")

test_that("Exact test and permutation test are consistent",{
  set.seed(1)
  tar_feat1 <- create_feature_target(10, 390, 0, 600) 
  tar_feat2 <- create_feature_target(9, 391, 1, 599)
  tar_feat3 <- create_feature_target(8, 392, 2, 598)
  fast.results <- test_features(tar_feat1[, 1, drop=FALSE], 
                     cbind(tar_feat1[,2], tar_feat2[,2], tar_feat3[,2]))
  
  m <- 10000
  perm.results <- test_features(tar_feat1[, 1, drop=FALSE], 
                cbind(tar_feat1[,2], tar_feat2[,2], tar_feat3[,2]),
                times = m, quick = FALSE)
  alfa <- 0.1
  conf.intervals <- sapply(perm.results, 
                           function(x) {
                             x <- (x*m+0.5*qnorm(alfa)^2)/(m+qnorm(alfa)^2)
                             x+qnorm(c(alfa/2,1-alfa/2))*sqrt(x*(1-x)/(m+qnorm(alfa)^2))})
  
  for (i in 1:3) {
    expect_true(fast.results[i] < conf.intervals[2, i])
    expect_true(fast.results[i] > conf.intervals[1, i])
  } 
})
