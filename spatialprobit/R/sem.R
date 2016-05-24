# Estimate the spatial error model (SEM)
# y = xB + u
# u = pWu + e, e ~ N(0,sige*V); V=diag(v1,v2,...vn)
#
# Code based on sem_g function; see Spatial Econometrics Toolbox
#
# @param y
# @param X
# @param W spatial weight matrix
# @param ndraw number of MCMC iterations
# @param burn.in  number of MCMC burn-in to be discarded
# @param thinning MCMC thinning factor, defaults to 1
# @param m number of burn.in sampling in inner Gibbs sampling
# @param prior list of prior settings (a1, a2, c, T, nu, d0) for rho ~ Beta(a1,a2), beta ~ N(c, T),
#    1/sige ~ Gamma(nu,d0) prior on sigma, default: nu=0,d0=0 (diffuse prior)
#
# @param start  start value for the chain; list of three components rho and beta and sigma.
# @param m
# @param showProgress
sem <- function(y, X, W, ndraw=1000, burn.in=100, thinning=1,
  prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12,
  nu=0, d0=0, lflag = 0),
  start=list(rho=0.75, beta=rep(0, ncol(X)), sige=1),
  m=10, showProgress=FALSE){

  #start timer
  timet <- Sys.time()

  n  <- nrow( X )            # number of observations
  n1 <- nrow( X )
  n2 <- nrow( W )
  k <- ncol( X )             # number of of parameters/exogenous variables
  I_n <- sparseMatrix(i=1:n, j=1:n, x=1) # sparse identity matrix
  if (is.null(colnames(X))) colnames(X) <- paste("x",1:k,sep="")

  if( n1 != n2 && n1 != n ){
    stop("sem: wrong size of spatial weight matrix W")
  }
  # check if we have a constant term in X
  ind <- match( n, apply(X,2,sum))
  if( is.na(ind) ){
    cflag <- 0
    p     <- k
  }else if( ind == 1 ){
    cflag <- 1
    p     <- k - 1
  }else{
    stop("sem: intercept term must be in first column of the X-matrix")
  }

  # MCMC sampling of beta
  rho  <- start$rho          # start value of row
  beta <- start$beta         # start value of parameters, prior value, we could also sample from beta ~ N(c, T)
  sige <- start$sige         # start value for sigma_e

  # conjugate prior beta ~ N(c, T)
  # parametrize, default to diffuse prior, for beta, e.g. T <- diag(k) * 1e12
  c <- rep(0, k)             # prior distribution of beta ~ N(c, T) : c = 0
  if (is.matrix(prior$T) && ncol(prior$T) == k && isSymmetric(prior$T) && det(prior$T) > 0) {
    T <- prior$T               # prior distribution of beta ~ N(c, T) : T = I_n --> diffuse prior
  } else {
    T <- diag(k)*1e12
  }
  # prior for sige ~ IG(nu, d0)
  if (is.numeric(prior$nu)) {
    nu <- prior$nu
  } else {
    nu <- 0
  }
  if (is.numeric(prior$d0)) {
    d0 <- prior$d0
  } else {
    d0 <- 0
  }

  V = rep(1, n)    # V=diag(v1,v2,...vn)
  inn = rep(1, n)  # initial value for V
  vi = inn

  TI <- solve(T)           # T^{-1}
  TIc <- TI%*%c            # T^{-1}c
  S <- I_n - rho * W
  Wy <- W %*% y
  WX <- W %*% X

  rmin       <- -1   # use -1,1 rho interval as default
  rmax       <-  1
  ldetflag   <-  0   # default to 1999 Pace and Barry MC determinant approx
  tmp <- sar_lndet(ldetflag, W, rmin, rmax)
  detval <- tmp$detval

  # Some precalculated quantities for drawing rho
  # rho ~ Beta(a1, a2) prior
  a1         <-  1.0
  a2         <-  1.0
  if(is.numeric(prior$a1)) a1 <- prior$a1
  if(is.numeric(prior$a2)) a2 <- prior$a2

  u        <- runif(thinning * ndraw + burn.in)   # u ~ U(0, 1)
  nrho     <- 2001
  nmk      <- (n-k)/2
  detval1  <- detval[,1]  # SW: avoid multiple accesses to detval[,1]
  detval2  <- detval[,2]

  # matrix to store the beta + sige + rho parameters for each iteration/draw
  B <- matrix(NA, ndraw, k+2)
  colnames(B) <- c(colnames(X), "sige", "rho")

  cc <- 0.2         # initial tuning parameter for M-H sampling
  acc <- 0          # number of accepted samples
  acc_rate <- rep(NA, thinning * ndraw + burn.in)

  # progress bar
  if (showProgress) {
    pb <- txtProgressBar(min=0, max=(thinning * ndraw + burn.in), initial=0, style=3)
  }

  # names of non-constant parameters
  if(cflag == 0) {
    namesNonConstantParams <- colnames(X)
  } else {
    namesNonConstantParams <- colnames(X)[-1]
  }

  for (i in (1 - burn.in):(ndraw * thinning)) {
  
# TODO: heteroskedastic model; e ~ N(sige* V), diag(V) = diag(v1,v2,...vn) <> (1, 1, ..., 1)
  
# homoskedastic model; e ~ N(sige* V) = N(sige*I); diag(V) = diag(v1,v2,...vn) = (1, 1, ..., 1)

  # update beta
  SX <- X - rho*WX
  AI <- solve(t(SX) %*% SX + sige * TI)
  Sy <- y - rho*Wy
  b <- t(SX) %*% Sy + sige * TIc
  b0 <- AI %*% b
  beta <- as.double(rmvnorm(n=1, mean=b0, sigma=as.matrix(sige*AI)))

  # update sige
  nu1 <- n + 2*nu
  e <- Sy - SX %*% beta
  d1 <- 2*d0 + crossprod(e)
  chi <- rchisq(n=1,df=nu1)
  sige <- as.double(d1/chi)
  
  # update rho using numerical integration
  #rho <- draw_rho(detval,y,x,Wy,Wx,V,n,k,rmin,rmax,rho);

  # update rho using metropolis-hastings
  xb <- X%*%beta;
  rhox <- c_sem(rho,y,X,beta,sige,I_n,W,detval1,detval2,rep(1,n),a1,a2);
  accept <- 0
  rho2 <- rho + cc * rnorm(1)
  while(accept == 0) {
    if ((rho2 > rmin) & (rho2 < rmax)) {
      accept <- 1
    }
    else {
      rho2 <- rho + cc * rnorm(1)
    }
  }
  rhoy <- c_sem(rho2,y,X,beta,sige,I_n,W,detval1,detval2,rep(1,n),a1,a2)
  ru <- runif(1,0,1) # TODO: Precalculate ru
  if ((rhoy - rhox) > exp(1)) {
    p <- 1
  } else {
    ratio <- exp(rhoy-rhox)
    p <- min(1,ratio)
  }
  if (ru < p) {
    rho <- rho2
    acc <- acc + 1
  }
  acc_rate[i + burn.in] <- acc/(i + burn.in)
  # update cc based on std of rho draws
  if (acc_rate[i + burn.in] < 0.4) {
    cc <- cc/1.1;
  }
  if (acc_rate[i + burn.in] > 0.6) {
    cc <- cc*1.1;
  }


  ##############################################################################

  # update S and H
  S <- I_n - rho * W
  H <- (1/sige) * t(S) %*% S      # H = t(S)%*%S  / SW: crossprod(S) does not seem to work!

  if (i > 0) {
    if (thinning == 1) {
      ind <- i
    }
    else if (i%%thinning == 0) {
      ind <- i%/%thinning
    } else {
      next
    }
    B[ind,] <- c(beta, sige, rho)

  ##############################################################################

  }

  if (showProgress) setTxtProgressBar(pb, i + burn.in) # update progress bar
  }

  if (showProgress)  close(pb) #close progress bar

  # fitted values for estimates (based on z rather than binary y like in fitted(glm.fit))
  # (on reponse scale y vs. linear predictor scale z...)
  beta  <- colMeans(B)[1:k]
  sige  <- colMeans(B)[k+1]
  rho   <- colMeans(B)[k+2]
  S     <- (I_n - rho * W)
  fitted.values   <- X %*% beta                     # E[z | beta] = (X * beta)
  fitted.response <- as.numeric(fitted.values >= 0)
  # TODO: linear.predictors  vs. fitted.values

  # result
  results       <- NULL
  results$time  <- Sys.time() - timet
  results$nobs  <- n          # number of observations
  results$nvar  <- k          # number of explanatory variables
  results$y     <- y
  results$zip   <- n - sum(y) # number of zero values in the y-vector
  results$beta  <- colMeans(B)[1:k]
  results$sige  <- sige
  results$rho   <- colMeans(B)[k+2]
  results$coefficients <- colMeans(B)
  results$fitted.values <- fitted.values
  #results$fitted.reponse <- fitted.reponse  # fitted values on reponse scale (binary y variable)
  results$ndraw <- ndraw
  results$nomit <- burn.in
  results$a1        <- a1
  results$a2        <- a2
  results$nu        <- nu
  results$d0        <- d0
  results$rmax      <- rmax
  results$rmin      <- rmin
  results$tflag     <- "plevel"
  results$lflag     <- ldetflag
  results$cflag     <- cflag
  results$lndet     <- detval
  results$names     <- c(colnames(X), "sige", "rho")
  results$B         <- B        # (beta, sige, rho) draws
  results$bdraw     <- B[,1:k]  # beta draws
  results$sdraw     <- B[,k+1]  # sige draws
  results$pdraw     <- B[,k+2]  # rho draws
  results$W <- W
  results$X <- X

  #results$predicted <- # prediction required. The default is on the scale of the linear predictors
  class(results)    <- "sem"
  return(results)
}