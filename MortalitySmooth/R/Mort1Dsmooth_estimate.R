Mort1Dsmooth_estimate <-
function(x, y, offset, wei,
                                  psi2, B, lambda,
                                  DtD, a.init,
                                  MON, TOL1, MAX.IT){
  ## Input:
  ## x: abcissae of data
  ## y: count response
  ## offset: an a priori known component
  ## wei: weigths
  ## psi2: overdispersion parameter
  ## B: B-splines basis
  ## lambda: smoothing parameter
  ## DtD: matrix of differences
  ## a.init: starting coefficients
  ## MON: logical on monitoring
  ## TOL1: relative convergence tolerance
  ## MAX.IT: the maximum number of iterations
  
  ## Output: a list containing
  ## a: fitted coefficients
  ## h: diag of the hat-matrix
  ## df: effective dimension
  ## aic: Akaike Information Criterion
  ## bic: Bayesian Information Criterion
  ## dev: Poisson-deviance
  ## tol: tolerance level
  ## BtWB: inner product of B-splines and weights
  ## P: penalty term

  ## penalty stuff
  lambdaP <- lambda * DtD
  ## initialize
  tol <- 1
  i <- 0
  a <- a.init
  a.old <- 10
  ## monitoring?
  if(MON){
    cat("lambda =", lambda, "\n")
    cat("Iter         tol", "\n")}
  while(tol > TOL1 && i < MAX.IT){
    i <- i+1
    ## update the coefficients
    a <- Mort1Dsmooth_update(x=x, y=y,
                             offset=offset, wei=wei,
                             psi2=psi2, B=B,
                             lambdaP=lambdaP, a=a)
    ## conputing the current tolerance level
    tol <- max(abs(a - a.old)/abs(a))
    ## replace the old coeff
    a.old <- a
    ## monitoring?
    if(MON){
      cat(i, "      ", tol, "\n")}
  }
  if(i > (MAX.IT-1)) {
    warning(paste("parameter estimates did NOT converge in",
                  MAX.IT,
            "iterations. Increase MAX.IT in control."))
  }
  ## final step after convergence
  eta <- B%*%a
  mu  <- exp(eta + offset)
  w   <- c(wei*mu)
  z   <- c(wei*((y - mu)/mu + eta))
  ## regression
  BtWB   <- t(B) %*% (w * B)
  P      <- psi2*lambdaP
  BtWB.P <- BtWB + P
  BtWz   <- t(B) %*% (w * z)
  a      <- solve(BtWB.P, BtWz)
  ## coefficients
  a <- matrix(a, ncol = 1)
  ## diag of the hat-matrix
  H <- solve(BtWB.P, BtWB)
  h <- diag(H)
  ## diagnostics
  ## replace zeros in response
  y1       <- y
  y1[y==0] <- 10^(-4)
  ## deviance
  dev <- 2*(sum(wei*(y1 * log(y1/mu) ), na.rm = TRUE))
  ## effective dimension
  df <- sum(h)
  ## Akaike Information Criterion
  aic <- dev/psi2 + 2 * df
  ## Bayesian Information Criterion
  bic <- dev/psi2 + log(sum(wei)) * df
  ## output
  llist <- list(a=a, h=h, 
                df=df, aic=aic, bic=bic,
                dev=dev, tol=tol,
                BtWB=BtWB, P=P)
  llist
}
