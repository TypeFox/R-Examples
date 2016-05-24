context("Testing consistency and timings of conjugate gradient algorithm - VERY PRELIMINARY")

test_that("dev_conjugate_gradient", {

  require(quadrupen)

  get.coef <- function(x,y) {
    lambda1 <- .25

    enet.ref <- elastic.net(x,y,lambda1=lambda1, control=list(timer=TRUE))
    enet.ref.bot <- elastic.net(x,y,lambda1=lambda1*2)

    enet.cg  <- elastic.net(x,y,lambda1=lambda1, control=list(timer=TRUE,usechol=FALSE,threshold=1e-3))
    enet.cg.warm <- elastic.net(x,y,lambda1=lambda1, beta0 = enet.ref.bot@coefficients, control=list(timer=TRUE,usechol=FALSE,threshold=1e-3))

    cat("\n\tTimings with warm-restart along the path")
    cat("\n\t\tfrom stratch (cholesky): ",enet.ref@monitoring$internal.timer)
    cat("\n\t\tfrom stratch (conjugate-gradient): ",enet.cg@monitoring$internal.timer)
    cat("\n\t\tCG starting from sparser solution: ",enet.cg.warm@monitoring$internal.timer)

    return(list(
      coef.ref=as.matrix(enet.ref@coefficients),
      coef.cg.warm=as.matrix(enet.cg.warm@coefficients),
      coef.cg =as.matrix(enet.cg@coefficients)))
  }

  ## PROSTATE DATA SET
  load("prostate.rda")
  x <- as.matrix(x)

  ## Run the tests...
  cat("\n  * tiny-size problem...")
  out <-get.coef(x,y)
  expect_equal(out$coef.cg.warm,out$coef.ref,tolerance=1e-3)
  expect_equal(out$coef.cg ,out$coef.ref    ,tolerance=1e-3)

  ## RANDOM DATA
  seed <- sample(1:10000,1)
  ## cat(" #seed=",seed)
  set.seed(seed)

  beta <- rep(rep(c(0,1,0,-1,0), c(25,10,25,10,25)),5)
  n <- 300
  p <- length(beta)

  mu <- 3 # intercept
  sigma <- 30 # huge noise
  Sigma <- matrix(0.95,p,p) # huge correlation
  diag(Sigma) <- 1

  x <- as.matrix(matrix(rnorm(p*n),n,p) %*% chol(Sigma))
  y <- 10 + x %*% beta + rnorm(n,0,10)

  ## Run the tests...
  cat("\n  * small-size problem, with correlation...")
  out <-get.coef(x,y)
  expect_equal(out$coef.cg.warm,out$coef.ref,tolerance=1e-3)
  expect_equal(out$coef.cg     ,out$coef.ref,tolerance=1e-3)

})

