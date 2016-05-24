context("W.combine")


## try and feed wrong things

test_that("W.combine throws errors with wrong param types",{
  expect_that(W.combine(1,1),throws_error(regexp="W must"))
  expect_that(W.combine(c(1,2,3),'foobar'),throws_error(regexp="W must"))
})

test_that("W.combine throws an error when param W and n are not equal length",{
  expect_that(W.combine(list(),1),throws_error(regexp="need equal"))
})

test_that("W.combine computes correctly",{
  ## this 
  W.t1<-list(rep(1,3),rep(2,3),rep(3,3))
  ## and this
  W.t2<-list(1,2,3)
  n<-c(1,1,1)
  ##should be equivalent
  expect_equal(unique(W.combine(W.t1,n)),W.combine(W.t2,n))
})

context("Z.value")

##try and feed wrong things

test_that("Z.value throws errors with wrong param types",{
  expect_that(Z.value(W=list(1,2),Wstar=c(1),n.in=1,n.out=1),throws_error(regexp="W, Wstar must"))
  expect_that(Z.value(W=list(1,2),Wstar=list(1,2),n.in=c(1),n.out=c(1,2)),throws_error(regexp="W, Wstar must"))
  ##currently no parameter checking if called in list context
})

test_that("Z.value function correctly in list and vector context",{
  ##test vector
  Z.vect<-Z.value(W=1,Wstar=c(1,2),n.in=1,n.out=2)
  expect_is(Z.vect,"list")
  expect_equal(names(Z.vect),c('Z.theoretical','Z.empirical'))
  expect_equal(Z.vect$Z.empirical$data,'1')
  ## test list
  Z.list<-Z.value(W=list(1,1),Wstar=list(c(1,2),c(1,2)),n.in=c(1,1),n.out=c(2,2))
  expect_equal(Z.list$Z.empirical$data,'0.5')
  ## vector and list are chosen so that they return same Z statistic
  expect_equal(Z.vect$Z.empirical$statistic,Z.list$Z.empirical$statistic)
  expect_equal(Z.vect$Z.theoretical$statistic,Z.list$Z.theoretical$statistic)
})

