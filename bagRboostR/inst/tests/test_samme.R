context("SAMME as described by Zhu")
test_that("samme returns a character vector of predicted classes",{
  outcomes <- c("laying", "sitting", "standing", "walking", "walking_downstairs", "walking_upstairs")
  train <- data.frame(activity=sample(outcomes,10,replace=T),V1=rnorm(10),V2=rnorm(10),V3=rnorm(10))
  test <- data.frame(V1=rnorm(10),V2=rnorm(10),V3=rnorm(10))
  p<-samme(activity~.,train,test,3,ntree=5,trace=F)
  expect_is(p,"character")
  expect_equal(length(p),nrow(test))
})

test_that("error returns the sum of the misclassified weights", {
  w <- rep(0.1,10)
  misses <- c(rep(1,2),rep(0,8)) 
  expect_equal(error(w,misses),0.2)
  w<-c(0.05,0.01,0.04,0.2,0.015,0.085,0.2,0.1,0.18,0.12)
  expect_equal(error(w,misses),0.06)
  misses <- c(1,0,1,0,1,0,1,0,1,0)
  expect_equal(error(w,misses),0.485)
})

test_that("alpha returns the log of (1-err)/err plus the log of one less the number of classes if err > 0 and <= 0.5", {
  err <- 0.05
  K <- 4
  expect_equal(round(alpha(err,K),6),4.043051)
})
  
test_that("alpha returns 1 if err > 0.5", {
  err <- 0.51
  expect_equal(alpha(err,K),1)
})

test_that("alpha returns 20 if err <= 0",{
  err <- 0
  expect_equal(alpha(err,K),20)  
})

test_that("weights returns the renormalized exp of the current alpha times misses multiplied by the previous weights", {
  w<-rep(0.1,10)
  misses <- c(rep(1,2),rep(0,8))
  err <- error(w,misses)
  K <- 6
  a <- alpha(err,K)
  expect_equal(round(weights(w,a,misses),8),c(rep(0.41666667,2),rep(0.02083333,8)))
})

test_that("renorm returns the normalized vector of weights", {
  tmp_w <- c(rep(5.459815,2),rep(0.1,8))
  expect_equal(renorm(tmp_w),c(rep(0.465869230,2),rep(0.008532693,8)))
})

test_that("the sum of weights after renorm is 1", {
  tmp_w <- c(rep(5.459815,2),rep(0.1,8))
  expect_equal(sum(renorm(tmp_w)),1)
})

test_that("finalPrediction returns a vector of the predicted classes", {
  C <- matrix(nrow=5,ncol=2)
  C[,1] <- c("standing","standing","sitting","sitting","sitting")
  C[,2] <- rep("standing",5)
  A <- c(6.535043, 6.527759)
  expect_equal(finalPrediction(C,A),c("standing","standing","sitting","sitting","sitting"))
  A <- c(6.527759,6.535043)
  expect_equal(finalPrediction(C,A),rep("standing",5))
})

test_that("samplePrediction returns the class with the maximum votes, weighted by classifier vote, ties broken at random", {
  sample <- c("sitting", "standing", "sitting", "sitting", "sitting")
  A <- c(6.535043, 6.527759, 6.904823, 7.690789, 6.657183)
  expect_equal(samplePrediction(sample,A),"sitting")
})

test_that("votesTable returns a table with outcomes and the sum of weighted votes",{
  sample <- c("sitting", "standing", "sitting", "sitting", "sitting")
  A <- c(6.535043, 6.527759, 6.904823, 7.690789, 6.657183)
  table <- unlist(list(sitting=27.78784,standing=6.52776))
  expect_equal(round(votesTable(sample,A),5),table)
})

test_that("weightedClassVote returns the sum of votes cast for a class k each multiplied by their classifier weight a", {
  k <- "standing"
  sample <- c("standing","standing","sitting","standing","standing")
  A <- c(6.709469, 6.659018, 6.855231, 7.323083, 6.259542)
  expect_equal(round(weightedClassVote(k,sample,A),5),26.95111)
  k <- "sitting"
  expect_equal(round(weightedClassVote(k,sample,A),6),6.855231)
})

test_that("prediction returns the class with the most classifier votes, tie broken at random", {
  votes.table <- unlist(list(standing=26.775759,sitting=7.529446))
  expect_equal(prediction(votes.table),"standing")
  votes.table <- unlist(list(standing=26.775759,sitting=7.529446,laying=26.775759))
  expect_true(prediction(votes.table) %in% c("standing","laying"))
})

test_that("maxVotes returns the class or classes that received the most classifier votes", {
  votes.table <- unlist(list(standing=26.775759,sitting=7.529446))
  expect_equal(maxVotes(votes.table),"standing")
  votes.table <- unlist(list(standing=26.775759,sitting=7.529446,laying=26.775759))
  expect_equal(maxVotes(votes.table),c("standing","laying"))
  votes.table <- unlist(list(standing=7.529446,sitting=7.529446,laying=26.775759))
  expect_equal(maxVotes(votes.table),"laying")
})