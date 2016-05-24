context("Consistency of penscale (vs theoretical and 'glmnet')")

test_that("weighted_ quad2theo", {

  ## SIMPLE CHECK: used identity for design matrix
  n <- 100
  p <- 100
  x <- diag(rep(1,n))
  y <- rnorm(100)

  ## no penscale...
  lasso.quad <- elastic.net(x, y, intercept=FALSE, lambda2=0)
  theo.path  <- t(sapply(lasso.quad@lambda1, function(lambda) y*pmax(0,1-lambda/abs(y))))
  expect_that(as.matrix(lasso.quad@coefficients), is_equivalent_to(theo.path))

  ## glmnet is not equal to what is expected... probably due to the intercept treatment
  ##  lasso.glmn <- glmnet(x,y, intercept=FALSE,lambda.min.ratio=1e-2, thresh=1e-20)
  ##  theo.path  <- t(sapply(lasso.glmn$lambda*sqrt(n), function(lambda) y*pmax(0,1-lambda/abs(y))))
  ##  expect_that(as.matrix(t(lasso.glmn$beta)), is_equivalent_to(theo.path))

  ## random penscale...
  w <- 1/runif(p,0.5,1)
  w <- w/sum(w)*p ## to fit glmnet rescaling
  lasso.quad <- elastic.net(x, y,intercept=FALSE, penscale = w, lambda2=0)
  theo.path  <- t(sapply(lasso.quad@lambda1, function(lambda) y*pmax(0,1-lambda*w/abs(y))))
  expect_that(as.matrix(lasso.quad@coefficients), is_equivalent_to(theo.path))

  ## glmnet with intercept fit with quadrupen and the theory
  w <- 1/runif(p,0.5,1)
  w <- w/sum(w)*p ## to fit glmnet rescaling
  lasso.glmn <- glmnet(x,y, penalty.factor=w,lambda.min.ratio=1e-2, thresh=1e-20)
  lasso.quad <- elastic.net(x,y, lambda1=lasso.glmn$lambda*sqrt(n), penscale = w, lambda2=0)
  expect_that(as.matrix(t(lasso.glmn$beta)), is_equivalent_to(as.matrix(lasso.quad@coefficients)))
  ## Check the intercept term also
  expect_that(lasso.glmn$a0, is_equivalent_to(lasso.quad@mu))

})

