########## R-function: aspGetTheta ##########

# Estimates the spline coefficients 
# of the variance of random effects

# Last changed: 30 JUL 2007


"aspGetTheta" <-  function (theta.info, niter.var, family, weights,spar.method) 
{

#get the model matrices and current model estimates

  Xb <- theta.info$model.matrices$Xb
  Zb <- theta.info$model.matrices$Zb
  Wc <- theta.info$model.matrices$Wc
  kc <- theta.info$model.matrices$kc
  kb <- theta.info$model.matrices$kb
  ZVZ <- theta.info$model.matrices$ZVZ
  XVX <- theta.info$model.matrices$XVX
  ZVX <- theta.info$model.matrices$ZVX

  sigma.theta <- theta.info$sigma.theta
  theta <- theta.info$theta
  sigma.eps <- theta.info$asp.info$fit$sigma^2
  sigma.b <- as.vector(sigma.eps * exp(2 * unlist(theta.info$asp.info$fit$modelStruct$reStruct)))
  sigma.b <- rep(sigma.b, kb)
  Y <- theta.info$asp.info$info$y
  beta <- theta.info$asp.info$fit$coef$fixed
  ind <- rep(c(TRUE, FALSE), length(sigma.theta))
  sigma.theta.d <- rep(1/sigma.theta, each = 2)
  sigma.theta.d[ind == TRUE] <- 0
  sigma.theta.d <- rep(sigma.theta.d, kc)
  D <- diag(sigma.theta.d)
  U <- Y
  resid <- t(U - Xb %*% beta)
  
#redefine matrices in case of non-normal response

  if (family != "gaussian")
    {
    mu.hat <- theta.info$asp.info$fit$fitted
    if (family == "poisson") 
      V <- c(mu.hat)
    if (family == "binomial") 
      V <- c(mu.hat * (1 - mu.hat) * weights)
    eta <- link(mu.hat, family)
    V12 <- sqrt(V)
    U <- c(eta + (Y - mu.hat)/V)
    Zb <- t(t(Zb) * V12)
    resid <- t((U - c(Xb %*% beta))*V12)
    Xb<- t(t(Xb)*V12)
    sigma.eps <- 1
  }

#iterate for the spline coeficients of the variance of random effects

  for (k in 1:niter.var) {
    Sigma.b.inverse <- exp(-c(Wc %*% theta))
    ridge <- chol(ZVZ + diag(Sigma.b.inverse * sigma.eps/sigma.b))
    Ridge.inv <- backsolve(ridge, diag(rep(1, nrow(ridge))))
    ridge.inv <- Ridge.inv %*% t(Ridge.inv)
    zz <- Zb %*% ridge.inv
    alpha <- c(resid %*% zz)
    ZZ <- t(Zb) %*% zz
    if (spar.method=="REML")
      {
        zzXb <- t(zz)%*%Xb
        PP <- (Sigma.b.inverse*sigma.eps/sigma.b)*(zzXb%*%solve(XVX-t(zzXb)%*%ZVX)%*%t(zzXb))
        ZZ <- ZZ-PP
     }
    wdf <- diag(ZZ)
    vdf <- diag(diag(ZZ %*% ZZ))
    wvw <- t(Wc) %*% vdf %*% Wc/2
    score <- diag(D)*theta-t(Wc) %*% (Sigma.b.inverse * (alpha^2)/sigma.b - wdf)/2
    ridge.theta <- chol(wvw + D)
    Ridge.theta.inv <-  backsolve(ridge.theta, diag(rep(1, nrow(ridge.theta))))
    fisher <-  Ridge.theta.inv%*%t(Ridge.theta.inv)
    theta1 <- theta-fisher %*% score
    epsilon.theta <- sum((theta - theta1)^2)/sum(theta^2)
    theta <- theta1
    if (epsilon.theta <= 1e-06) 
      break
    if (k == niter.var) 
      stop("Iteration limit reached without convergence in variance parameter")
  }

#get the variance of the spline coeficients of the variance of random effects

  ind1 <- rep(rep(c(FALSE,TRUE),length(kc)/2),kc)
  zvz <- wvw[ind1,ind1]
  D <- diag(D)
  D[D != 0] <- 1
  c <- (D * theta)^2
  c <- c[c != 0]
  kcc <- kc[ind == FALSE]
  ind2 <- rep(1:length(kcc), kcc)
  ctc <- as.vector(by(c, ind2, sum))
  zvz.df <- solve(zvz+diag(rep(1/sigma.theta,kcc)))%*%zvz
  zvz.df <- as.vector(by(diag(zvz.df), ind2, sum))
  sigma.theta <- ctc/zvz.df
 
  return(list(theta = theta, sigma.theta = sigma.theta))
}
