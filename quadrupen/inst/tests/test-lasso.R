context("Consistency of the Lasso solution paths (package 'lars' and 'glmnet')")

test_that("lasso_quad2lars", {

  require(lars)

  get.lars <- function(x,y,intercept,normalize) {
      lasso.larsen <- lars(x,y,intercept=intercept,normalize=normalize)
      iols <- nrow(lasso.larsen$beta) ## remove last entry corresponding to the OLS estimator
      lambda1 <-  lasso.larsen$lambda ## usde the lars lambda grid
      lasso.quadru <- elastic.net(x,y, intercept=intercept, normalize=normalize,
                                  lambda1=lambda1, lambda2=0, control=list(method="quadra"))
      quad <- list(coef   = as.matrix(lasso.quadru@coefficients),
                   meanx  = lasso.quadru@meanx,
                   normx  = lasso.quadru@normx,
                   rss    = deviance(lasso.quadru))

      lars <- list(coef   = lasso.larsen$beta[-iols, ],
                   meanx  = lasso.larsen$meanx,
                   normx  = lasso.larsen$normx,
                   rss    = lasso.larsen$RSS[-iols])

      return(list(quad=quad,lars=lars))
  }

  ## PROSTATE DATA SET
  load("prostate.rda")
  x <- as.matrix(x)

  ## Run the tests...
  with.intercept <-get.lars(x,y,TRUE,TRUE)
  expect_that(with.intercept$quad,
              is_equivalent_to(with.intercept$lars))

  with.intercept.unnormalized <-get.lars(x,y,TRUE,FALSE)
  expect_that(with.intercept.unnormalized$quad,
              is_equivalent_to(with.intercept.unnormalized$lars))

  without.intercept <-get.lars(x,y,FALSE,TRUE)
  expect_that(without.intercept$quad,
              is_equivalent_to(without.intercept$lars))

  without.intercept.unnormalized <-get.lars(x,y,FALSE,FALSE)
  expect_that(without.intercept.unnormalized$quad,
              is_equivalent_to(without.intercept.unnormalized$lars))

  ## RANDOM DATA
  seed <- sample(1:10000,1)
  ## cat("\n#seed=",seed)
  set.seed(seed)

  beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
  n <- 100
  p <- length(beta)

  mu <- 3 # intercept
  sigma <- 30 # huge noise
  Sigma <- matrix(0.95,p,p) # huge correlation
  diag(Sigma) <- 1

  x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
  y <- 10 + x %*% beta + rnorm(n,0,10)

  ## Run the tests...
  with.intercept <-get.lars(x,y,TRUE,TRUE)
  expect_that(with.intercept$coef.quad,
              is_equivalent_to(with.intercept$coef.lars))

  with.intercept.unnormalized <-get.lars(x,y,TRUE,FALSE)
  expect_that(with.intercept.unnormalized$coef.quad,
              is_equivalent_to(with.intercept.unnormalized$coef.lars))

  without.intercept <-get.lars(x,y,FALSE,TRUE)
  expect_that(without.intercept$coef.quad,
              is_equivalent_to(without.intercept$coef.lars))

  without.intercept.unnormalized <-get.lars(x,y,FALSE,FALSE)
  expect_that(without.intercept.unnormalized$coef.quad,
              is_equivalent_to(without.intercept.unnormalized$coef.lars))

})

test_that("lasso_quad2glmnet", {

  require(glmnet)

  ## SECOND CHECK: compare to glmnet with prescaling of x
  x <- matrix(rnorm(100*50),100,50)
  y <- rnorm(100)
  y <- y-mean(y)
  n <- nrow(x)
  p <- ncol(x)

  ## If thresh is set to the default, the test won't pass!!!
  ## This is beacause coordinate descent is fast yet not extremely accurate
  lasso.glmn <- glmnet(x,y, lambda.min.ratio=1e-2, thresh=1e-20)
  lasso.quad <- elastic.net(x,y, lambda1=lasso.glmn$lambda*sqrt(n), lambda2=0)

  quad <- list(coef   = as.matrix(lasso.quad@coefficients),
               mu     = lasso.quad@mu,
               fitted = as.matrix(lasso.quad@fitted))

  glmn <- list(coef   = as.matrix(t(lasso.glmn$beta)),
               mu     = lasso.glmn$a0,
               fitted = predict(lasso.glmn,x))


  expect_that(quad, is_equivalent_to(glmn))

})
