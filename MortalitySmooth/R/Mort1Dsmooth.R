Mort1Dsmooth <-
function(x, y, offset, w,
                         overdispersion=FALSE,
                         ndx=floor(length(x)/5),
                         deg=3, pord=2, 
                         lambda=NULL, df=NULL, method=1,
                         coefstart=NULL,
                         control=list()){
  ## Input:
  ## x: abcissae of data
  ## y: count response
  ## offset: an a priori known component (optional)
  ## w: a vector of weights to be used in the
  ##    fitting process (optional)
  ## overdispersion: logical on the presence of an
  ##                 overdispersion parameter.
  ##                 Default: FALSE
  ## ndx: number of internal knots -1.
  ##      Default: floor(length(x)/5)
  ## deg: degree of the B-splines. Default: 3
  ## pord: order of differences. Default: 2
  ## lambda: smoothing parameter (optional)
  ## df: a number which specifies the degrees
  ##     of freedom (optional)
  ## method: the method for controlling
  ##         the amount of smoothing. Default: 1
  ## coefstart: eventual starting coefficients
  ## control: a list of control parameters
  ## MON: logical on monitoring
  ## TOL1: convergence of the IWLS algorithm.
  ##       Default: 1e-06.
  ## TOL2: difference between two adjacent
  ##       smoothing parameters in the grid search,
  ##       log-scale. Default: 0.1.
  ## RANGE: range in between lambda is searched, log-scale.
  ##        Default: [10^-8 ; 10^8]
  ## MAX.IT: the maximum number of iterations. Default: 50
  


  ## Output: a Mort1Dsmooth object containing
  ## aic: Akaike Information Criterion
  ## bic: Bayesian Information Criterion
  ## dev: Poisson-deviance
  ## lev: diag of the hat-matrix
  ## df: effective dimension
  ## lambda: the selected (given) lambda
  ## psi2: (estimated) overdispersion parameter
  ## ndx: number of internal knots -1
  ## deg: degree of the B-splines
  ## pord: order of differences
  ## B: B-splines basis
  ## x: abcissae of data
  ## y: count responses
  ## offset: an a priori known component
  ## w: (only for weighted fits) the specified weights
  ## coef: fitted (penalized) B-splines coefficients  
  ## residuals: the deviance residuals
  ## fitted.values: fitted counts
  ## linear.predictor: fitted linear predictor
  ## logmortality: fitted mortality rates in log-scale
  ## call: the matched call
  ## n: number of observations
  ## tolerance: convergence tolerance

  ## checkings:
  if(missing(w)){
    w <- rep(1, length(y))
  }
  ## about x and y:
  if(missing(y)){
    if(is.list(x)){
      if(any(is.na(match(c("x", "y"), names(x)))))
        stop("cannot find x and y in list")
      x <- x$x
      y <- x$y
    }
    else if(is.complex(x)) {
            y <- Im(x)
            x <- Re(x)
          }
    else if(is.matrix(x) && ncol(x) == 2) {
      y <- x[, 2]
      x <- x[, 1]
    }
    else {
      y <- x
      x <- time(x)
    }
  }
  ## about the offset
  m <- length(y)
  if(missing(offset)) {
    offset <- rep(0, m)
  }else{
    if(length(offset) == 1) {
      offset <- rep(offset, m)
    }
    else{
      offset=offset
    }
  }  
  check <- Mort1Dsmooth_checker(x=x, y=y,
                                offset=offset, w=w,
                                overdispersion=overdispersion,
                                ndx=ndx, deg=deg,
                                pord=pord, 
                                lambda=lambda, df=df,
                                method=method,
                                coefstart=coefstart,
                                control=control)
  x <- check$x
  y <- check$y
  m <- check$m
  offset <- check$offset
  offsetINIT <- check$offsetINIT
  wei <- check$w
  over <- check$overdispersion
  ndx <- check$ndx
  deg <- check$deg
  pord <- check$pord
  lambda <- check$lambda
  df <- check$df
  MET <- check$method
  a.init <- check$coefstart
  MON <- check$control$MON
  TOL1 <- check$control$TOL1
  TOL2 <- check$control$TOL2
  RANGE <- check$control$RANGE
  MAX.IT <- check$control$MAX.IT
  call <- match.call()
  ## B-splines basis
  xl <- min(x)
  xr <- max(x)
  xmax <- xr + 0.01 * (xr - xl)
  xmin <- xl - 0.01 * (xr - xl)
  B <- MortSmooth_bbase(x, xmin, xmax, ndx, deg)
  ## penalty stuff
  nb <- ncol(B)
  D. <- diff(diag(nb), diff=pord)
  DtD <- t(D.)%*%D.
  ## General initialize:
  if(is.null(a.init)){
    y[is.na(y)] <- 0
    ## 0) simple poisson-GLM with ages(years)
    ##    only for the interpolation cases
    fit0 <- glm(round(y) ~ x + offset(offset),
                family=poisson, weights=wei)
    ## 1) simple penalized-poisson-GLM with B
    etaGLM <- log(fit0$fitted) - offset
    eta0 <- log((y + 1)) - offset
    eta0[wei==0] <- etaGLM[wei==0]
    mu0 <- exp(offset + eta0)
    w0 <- wei*mu0
    z0 <- wei*((y - mu0)/mu0 + eta0)
    BtWB <- t(B) %*% (w0 * B)
    BtWz <- t(B) %*% (w0 * z0)
    a.init <- solve(BtWB + 1e08 * DtD, BtWz)
  }else{
    a.init=a.init
  }
  psi2 <- 1

  ## ## plotting starting coef
  ## ran <- range(log(y/e), B%*%a.init,
  ##              na.rm=TRUE, finite = TRUE)
  ## plot(x, log(y/e), ylim=ran)
  ## lines(x, B%*%a.init, col=2, lwd=2)

  ## optimize AIC or BIC
  if(MET==1|MET==2){
    ## if overdisperion is true
    if(over){
      tol.over <- 10
      i.over <- 0
      while(tol.over > 1e-03 && i.over < 5){
        i.over <- i.over+1
        lambda.hat <- Mort1Dsmooth_optimize(x=x,
                                            y=y,
                                            offset=offset,
                                            wei=wei,
                                            psi2=psi2,
                                            B=B, DtD=DtD,
                                            a.init=a.init,
                                            MON=MON,
                                            TOL1=TOL1,
                                            TOL2=TOL2,
                                            RANGE=RANGE,
                                            MAX.IT=MAX.IT,
                                            MET=MET)
        FIT <- Mort1Dsmooth_estimate(x=x, y=y,
                                     offset=offset,
                                     wei=wei,
                                     psi2=psi2,
                                     B=B, 
                                     lambda=lambda.hat,
                                     DtD=DtD,
                                     a.init=a.init, 
                                     MON=MON, TOL1=TOL1,
                                     MAX.IT=MAX.IT)
        ## recalculating overdispersion parameter
        psi2.old <- psi2
        psi2 <- FIT$dev / (m - FIT$df)
        tol.over <- abs(psi2 - psi2.old)/abs(psi2)
      }
    }else{## if psi2==1
      lambda.hat <- Mort1Dsmooth_optimize(x=x,
                                          y=y,
                                          offset=offset,
                                          wei=wei,
                                          psi2=psi2,
                                          B=B, DtD=DtD,
                                          a.init=a.init,
                                          MON=MON,
                                          TOL1=TOL1,
                                          TOL2=TOL2,
                                          RANGE=RANGE,
                                          MAX.IT=MAX.IT,
                                          MET=MET)
      FIT <- Mort1Dsmooth_estimate(x=x, y=y, offset=offset,
                                   wei=wei,
                                   psi2=psi2,
                                   B=B, 
                                   lambda=lambda.hat,
                                   DtD=DtD,
                                   a.init=a.init, 
                                   MON=MON, TOL1=TOL1,
                                   MAX.IT=MAX.IT)
      psi2 <- FIT$dev / (m - FIT$df)
    }
    if(log10(lambda.hat)>=log10(RANGE[2]) |
       log10(lambda.hat)<=log10(RANGE[1])) {
 warning(paste("optimal lambda at the edge of the grid."))
    }
  }
  ## given lambda
  if(MET==3){
    lambda.hat <- lambda
    FIT <- Mort1Dsmooth_estimate(x=x, y=y, offset=offset,
                                 wei=wei,
                                 psi2=psi2,
                                 B=B, 
                                 lambda=lambda.hat, DtD=DtD,
                                 a.init=a.init, 
                                 MON=MON, TOL1=TOL1,
                                 MAX.IT=MAX.IT)
    psi2 <- FIT$dev / (m - FIT$df)
  }
  ## optimize given df
  if(MET==4){
    Mort1Dsmooth_opt_df <- function(X){
      FIT <- Mort1Dsmooth_estimate(x=x, y=y, offset=offset,
                                   wei=wei,
                                   psi2=psi2,
                                   B=B, 
                                   lambda=X, DtD=DtD,
                                   a.init=a.init, 
                                   MON=MON, TOL1=TOL1,
                                   MAX.IT=MAX.IT)
      return(abs(FIT$df - df))
    }
    by.lambda <- length(seq(log10(RANGE[1]),
                            log10(RANGE[2]),by=TOL2))
    lambda.hat <- cleversearch(fn=Mort1Dsmooth_opt_df,
                               lower=log10(RANGE[1]),
                               upper=log10(RANGE[2]),
                               ngrid=by.lambda,
                               startvalue=1,
                               logscale=TRUE,
                               verbose=FALSE)[[1]]
    if(log10(lambda.hat)>=log10(RANGE[2]) |
       log10(lambda.hat)<=log10(RANGE[1])){
  warning(paste("optimal lambda at the edge of the grid."))
    }
    FIT <- Mort1Dsmooth_estimate(x=x, y=y, offset=offset,
                                 wei=wei,
                                 psi2=psi2,
                                 B=B, 
                                 lambda=lambda.hat,
                                 DtD=DtD,
                                 a.init=a.init, 
                                 MON=MON, TOL1=TOL1,
                                 MAX.IT=MAX.IT)
    psi2 <- FIT$dev / (m - FIT$df)
  }
  aic <- FIT$aic
  bic <- FIT$bic
  df <- FIT$df
  dev <- FIT$dev
  coef <- FIT$a
  psi2 <- psi2
  h <- FIT$h
  tolerance <- FIT$tol
  eta.hat <- B%*%coef
  logmortality <- eta.hat
  fitted.values <- exp(eta.hat + offset)
  res <- sign(y - fitted.values) *
         sqrt(2 * (y * log(ifelse(y == 0, 1,
                                  y/fitted.values)) -
                   (y - fitted.values)))
  lin <- as.vector(eta.hat) + as.vector(offset)
  ## output
  object <- list(## fitted values
                 coefficients=as.vector(coef),
                 residuals=as.vector(res),
                 fitted.values=as.vector(fitted.values),
                 linear.predictors=lin,
                 logmortality=logmortality,
                 ## diagnostics
                 lev=h, df=df, deviance=dev,
                 aic=aic, bic=bic,
                 psi2=psi2, lambda=lambda.hat,
                 ## general
                 call=call, n=m, tolerance=tolerance,
                 ## smoothing specifications
                 ndx=ndx, deg=deg, pord=pord, B=B,
                 ## overdispersion=over,
                 ## observed values
                 x=x, y=as.vector(y),
                 offset=as.vector(offsetINIT),
                 w=as.vector(wei)
                 )
  class(object) <- "Mort1Dsmooth"
  object
}
