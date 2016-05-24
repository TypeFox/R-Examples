Mort2Dsmooth_estimate <-
function(x, y, Z, offset,
                                  psi2, wei, Bx, By,
                                  nbx, nby, RTBx, RTBy,
                                  lambdas, Px,
                                  Py, a.init,
                                  MON, TOL1, MAX.IT){
  ## Input:
  ## x: abcissae of data
  ## y: ordinate of data
  ## Z: matrix of count response
  ## offset: a matrix of an a priori known component
  ## psi2: overdispersion parameter
  ## wei: a matrix of weights
  ## Bx: B-splines basis for the x-axis
  ## By: B-splines basis for the y-axis
  ## nbx: number of B-splines for the x-axis
  ## nby: number of B-splines for the y-axis
  ## RTBx: tensors product of B-splines basis for the x-axis
  ## RTBy: tensors product of B-splines basis for the y-axis
  ## lambdas: smoothing parameters
  ## Px: penalty term/matrix for the x-axis
  ## Py: penalty term/matrix for the y-axis
  ## a.init: coefficients (in a matrix)
  ## MON: logical on monitoring
  ## TOL1: vector of the two relative convergence
  ##       tolerances for the smoothing parameters
  ## MAX.IT: the maximum number of iterations
  
  ## Output: a list containing
  ## a: fitted coefficients
  ## h: diag of the hat-matrix
  ## df: effective dimension
  ## aic: Akaike Information Criterion
  ## bic: Bayesian Information Criterion
  ## dev: Poisson-deviance
  ## tol: tolerance level
  
  ## penalty stuff
  P <- (lambdas[1] * Px) + (lambdas[2] * Py)
  ## initialize
  tol <- 1
  i <- 0
  a <- a.init
  a.old <- 10
  ## monitoring?
  if(MON){
    cat("lambda.x =", lambdas[1], "\n")
    cat("lambda.y =", lambdas[2], "\n")
    cat("Iter         tol", "\n")}
  while(tol > TOL1 && i < MAX.IT){
    i <- i+1
    ## update the coefficients
    a <- Mort2Dsmooth_update(x=x, y=y, Z=Z, offset=offset,
                             psi2=psi2, wei=wei, 
                             Bx=Bx, By=By,
                             nbx=nbx, nby=nby,
                             RTBx=RTBx, RTBy=RTBy, P=P, a=a)
    ## conputing the current tolerance level
    tol <- max(abs(a - a.old)/abs(a))
    ## replace the old coeff
    a.old <- a
    ## monitoring?
    if(MON){
      cat(i, "      ", tol, "\n")}
  }
  if(i > (MAX.IT-1)) {
    warning(paste("parameter estimates did NOT converge in", MAX.IT, "iterations. Increase MAX.IT in control."))
  }
  ## final step after convergence
  eta <- MortSmooth_BcoefB(Bx, By, a)
  mu <- exp(offset + eta)
  W <- mu
  z <- eta + (1/mu)* (Z - mu)
  z[which(wei==0)] <- 0
  WW <- wei*W
  BWB <- MortSmooth_BWB(RTBx, RTBy, nbx, nby, WW)
  BWBpP <- BWB + psi2*P
  BWz <- MortSmooth_BcoefB(t(Bx),t(By),(WW*z))
  a0 <- solve(BWBpP, c(BWz))
  a <- matrix(a0, nrow = nbx)
  ## diag of the hat-matrix
  H <- solve(BWBpP, BWB)
  h <- diag(H)
  ## diagnostics
  ## replace zeros in response
  Z1 <- Z
  Z1[Z==0] <- 10^(-4)
  ## deviance
  dev <- 2*(sum(wei*(Z1 * log(Z1/mu) ), na.rm = TRUE))
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
                BWB=BWB, P=P)
  llist
}
