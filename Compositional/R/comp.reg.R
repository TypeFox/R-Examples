################################
#### Regression for compositional data based on the log-ratio transformation
#### Tsagris Michail 6/2014
#### mtsagris@yahoo.gr
#### References: John Aitchison (2003)
#### The Statistical Analysis of Compositional Data p. 158-160 Blackburn Press
################################

comp.reg <- function(y, x, type = "classical", xnew = NULL, yb = NULL) {
  ## y is dependent variable, the compositional data
  ## x is the independent variable(s)
  ## type takes three values, either 'classical' or
  ## 'spatial' for spatial median regression.
  y <- as.matrix(y)
  y <- y/rowSums(y)  ## makes sure y is compositional data
  x <- as.matrix(x)

  ## alr transformation with the first component being the base
  if ( is.null(yb) )  {
    z <- log( y[, -1] / y[, 1] )
  } else {
    z <- yb
  }

  if (type == "classical") {
    runtime <- proc.time()
    mod <- multivreg(z, x, plot = FALSE, xnew = xnew)  ## classical multivariate regression
    res <- mod$suma
    di <- ncol(z)
    beta <- seb <- matrix(nrow = ncol(x) + 1, ncol = di)
    for (i in 1:di) {
     beta[, i] <- res[, 1, i]
     seb[, i] <- res[, 2, i]
    }
    rownames(seb) <- rownames(beta) <- rownames(res[, , 1])
    colnames(seb) <- colnames(beta) <- colnames(mod$fitted)
    est1 <- mod$est
    runtime <- proc.time() - runtime
  }

  if (type == "spatial") {
    mod <- spatmed.reg(z, x, xnew = xnew)  ## spatial median regression
    beta <- mod$beta
    seb <- mod$seb
    est1 <- mod$est
    runtime <- mod$runtime
  }

  est2 <- cbind(1, exp(est1))
  est <- est2/rowSums(est2)
  list(runtime = runtime, beta = beta, seb = seb, est = est)
}
