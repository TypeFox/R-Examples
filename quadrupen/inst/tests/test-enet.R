context("Consistency of the Elastic-net solution path (package 'elasticnet')")

test_that("enet_quad2elasticnet", {

  require(elasticnet)

  get.enet <- function(x,y,intercept,normalize=TRUE,naive=FALSE) {
    lambda2 <- runif(1,0,10)
    enet.larsen <- enet(x,y,lambda=lambda2,intercept=intercept,normalize=normalize)
    iols <- length(enet.larsen$penalty)
    lambda1 <- enet.larsen$penalty[-iols]/2
    enet.quadru <- elastic.net(x,y,intercept=intercept,normalize=normalize,
                               lambda1=lambda1, lambda2=lambda2, naive=naive)

    quad <- list(coef   = as.matrix(enet.quadru@coefficients),
                 meanx  = enet.quadru@meanx,
                 normx  = enet.quadru@normx)

    enet <- list(coef   = predict(enet.larsen, type="coefficients",naive=naive)$coefficients[-iols,],
                 meanx  = enet.larsen$meanx,
                 normx  = enet.larsen$normx)

    return(list(quad=quad, enet = enet))
  }

  ## PROSTATE DATA SET
  load("prostate.rda")
  x <- as.matrix(x)

  ## Run the tests...
  with.intercept <-get.enet(x,y,intercept=TRUE,naive=TRUE)
  expect_that(with.intercept$quad,
              is_equivalent_to(with.intercept$enet))

  without.intercept <- get.enet(x,y,intercept=FALSE,naive=TRUE)
  expect_that(without.intercept$quad,
              is_equivalent_to(without.intercept$enet))

  with.intercept <- get.enet(x,y,intercept=TRUE,naive=FALSE)
  expect_that(with.intercept$quad,
              is_equivalent_to(with.intercept$enet))

  without.intercept <-get.enet(x,y,intercept=FALSE,naive=FALSE)
  expect_that(without.intercept$quad,
              is_equivalent_to(without.intercept$enet))

  with.intercept <-get.enet(x,y,intercept=TRUE,normalize=FALSE,naive=TRUE)
  expect_that(with.intercept$quad,
              is_equivalent_to(with.intercept$enet))

  without.intercept <-get.enet(x,y,intercept=FALSE,normalize=FALSE,naive=TRUE)
  expect_that(without.intercept$quad,
              is_equivalent_to(without.intercept$enet))

  with.intercept <-get.enet(x,y,intercept=TRUE,normalize=FALSE,naive=FALSE)
  expect_that(with.intercept$quad,
              is_equivalent_to(with.intercept$enet))

  without.intercept <-get.enet(x,y,intercept=FALSE,normalize=FALSE,naive=FALSE)
  expect_that(without.intercept$quad,
              is_equivalent_to(without.intercept$enet))

  ## RANDOM DATA
  seed <- sample(1:10000,1)
  ## cat(" #seed=",seed)
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
  with.intercept <-get.enet(x,y,intercept=TRUE,naive=TRUE)
  expect_that(with.intercept$quad,
              is_equivalent_to(with.intercept$enet))

  without.intercept <-get.enet(x,y,intercept=FALSE,naive=TRUE)
  expect_that(without.intercept$quad,
              is_equivalent_to(without.intercept$enet))

  with.intercept <-get.enet(x,y,intercept=TRUE,naive=FALSE)
  expect_that(with.intercept$quad,
              is_equivalent_to(with.intercept$enet))

  without.intercept <-get.enet(x,y,intercept=FALSE,naive=FALSE)
  expect_that(without.intercept$quad,
              is_equivalent_to(without.intercept$enet))

  with.intercept <-get.enet(x,y,intercept=TRUE,normalize=FALSE,naive=TRUE)
  expect_that(with.intercept$quad,
              is_equivalent_to(with.intercept$enet))

  without.intercept <-get.enet(x,y,intercept=FALSE,normalize=FALSE,naive=TRUE)
  expect_that(without.intercept$quad,
              is_equivalent_to(without.intercept$enet))

  with.intercept <-get.enet(x,y,intercept=TRUE,normalize=FALSE,naive=FALSE)
  expect_that(with.intercept$quad,
              is_equivalent_to(with.intercept$enet))

  without.intercept <-get.enet(x,y,intercept=FALSE,normalize=FALSE,naive=FALSE)
  expect_that(without.intercept$quad,
              is_equivalent_to(without.intercept$enet))
  ## Run the tests...
  with.intercept <-get.enet(x,y,intercept=TRUE)
  expect_that(with.intercept$quad,
              is_equivalent_to(with.intercept$enet))

  without.intercept <-get.enet(x,y,intercept=FALSE)
  expect_that(without.intercept$quad,
              is_equivalent_to(without.intercept$enet))

})

