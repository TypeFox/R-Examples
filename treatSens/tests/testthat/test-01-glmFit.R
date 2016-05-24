context("glm fit")

test_that("C-glm implementation matches R's", {

  p <- 3
  n <- 50
  
  set.seed(0)
  offset <- rep(0.1, n)
  x <- matrix(c(rep(1, n), runif(n * (p - 1))), n)
  beta <- c(0.2, -0.4, 0.3)
  
  y <- rbinom(n, 1, pnorm(x %*% beta))
  
  x.1 <- x[,2]
  x.2 <- x[,3]
  
  glmCoefs <- as.numeric(glm(y ~ 1 + x.1 + x.2, offset = offset, family = binomial(link = probit))$coef)
  
  cibartCoefs <- .Call("treatSens_glmFit", y, NULL, x, NULL, offset)
  
  expect_equal(glmCoefs, cibartCoefs)
})
