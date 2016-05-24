################################
#### OLS regression for compositional data
#### Tsagris Michail 10/2014
#### mtsagris@yahoo.gr
#### References: Murteira, Jose M.R. and Ramalho, Joaquim J.S. (2013)
#### Regression analysis of multivariate fractional data
#### Econometric Reviews (to appear)
################################

ols.compreg <- function(y, x, B = 1, ncores = 1, xnew = NULL) {
  ## y is dependent variable, the compositional data
  ## x is the independent variable(s)
  ## B is the number of bootstrap samples used to obtain
  ## standard errors for the betas
  ## if B==1 no bootstrap is performed and no standard errors are reported
  ## if ncores=1, then 1 processor is used, otherwise
  ## more are used (parallel computing)

  runtime <- proc.time()
  y <- as.matrix(y)
  y <- y/rowSums(y)  ## makes sure y is compositional data
  x <- as.matrix( cbind(1, x) )
  d <- ncol(y) - 1  ## dimensionality of the simplex
  n <- nrow(y)  ## sample size
  z <- list(y = y, x = x)
  olsreg <- function(para, z) {
    y <- z$y
    x <- z$x
    d <- ncol(y) - 1
    be <- matrix(para, byrow = TRUE, ncol = d)
    mu1 <- cbind(1, exp(x %*% be))
    mu <- mu1 / rowSums(mu1)
    sum( (y - mu)^2 )
  }
  ## the next lines minimize the reg function and obtain the estimated betas
  ini <- as.vector( t( coef(lm(y[, -1] ~ x[, -1])) ) )  ## initial values
  options (warn = -1)
  qa <- nlm(olsreg, ini, z = z)
  qa <- nlm(olsreg, qa$estimate, z = z)
  qa <- nlm(olsreg, qa$estimate, z = z)
  beta <- matrix(qa$estimate, byrow = TRUE, ncol = d)
  seb <- NULL
  runtime <- proc.time() - runtime

  if (B > 1) {
  nc <- ncores
    if (nc == 1) {
      runtime <- proc.time()
      betaboot <- matrix(nrow = B, ncol = length(ini))
      for (i in 1:B) {
        ida <- sample(1:n, n, replace = TRUE)
        yb <- y[ida, ]
        xb <- x[ida, ]
        zb <- list(y = yb, x = xb)
        ini <- as.vector( t( coef(lm(yb[, -1] ~ xb[, -1])) ) )  ## initial values
        qa <- nlm(olsreg, ini, z = zb)
        qa <- nlm(olsreg, qa$estimate, z = zb)
        qa <- nlm(olsreg, qa$estimate, z = zb)
        betaboot[i, ] <- qa$estimate
      }
     s <- apply(betaboot, 2, sd)
     seb <- matrix(s, byrow = TRUE, ncol = d)
     runtime <- proc.time() - runtime

    } else {
      runtime <- proc.time()
      betaboot <- matrix(nrow = B, ncol = length(ini) )
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ww <- foreach::foreach(i = 1:B, .combine = rbind, .export="olsreg") %dopar% {
        ida <- sample(1:n, n, replace = TRUE)
        yb <- y[ida, ]
        xb <- x[ida, ]
        zb <- list(y = yb, x = xb)
        ini <- as.vector( t( coef(lm(yb[, -1] ~ xb[, -1])) ) )  ## initial values
        qa <- nlm(olsreg, ini, z = zb)
        qa <- nlm(olsreg, qa$estimate, z = zb)
        qa <- nlm(olsreg, qa$estimate, z = zb)
        betaboot[i, ] <- qa$estimate
      }
      stopCluster(cl)
      s <- apply(ww, 2, sd)
      seb <- matrix(s, byrow = TRUE, ncol = d)
      runtime <- proc.time() - runtime
    }
  }

  if ( is.null(xnew) ) {
    mu <- cbind( 1, exp(x %*% beta) )
    est <- mu / rowSums(mu)
  } else {
    xnew <- cbind(1, xnew)
    xnew <- as.matrix(xnew)
    mu <- cbind(1, exp(x %*% beta))
    est <- mu / rowSums(mu)
  }

  if ( is.null(colnames(x)) ) {
    p <- ncol(x) - 1
    rownames(beta) <- c("constant", paste("X", 1:p, sep = "") )
    if ( !is.null(seb) )  rownames(seb) <- c("constant", paste("X", 1:p, sep = "") )
  } else {
    rownames(beta)  <- c("constant", colnames(x)[-1] )
    if  ( !is.null(seb) ) rownames(seb) <- c("constant", colnames(x)[-1] )
  }

  list(runtime = runtime, beta = beta, seb = seb, est = est)
}
