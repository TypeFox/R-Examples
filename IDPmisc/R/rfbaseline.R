`rfbaseline` <-
function(x, y, span=2/3, NoXP=NULL, maxit=c(2,2), b=3.5,
		      weight=NULL, Scale=function(r) median(abs(r))/0.6745,
		      delta=NULL, SORT=TRUE, DOT=FALSE, init=NULL)
{
  if(!is.null(NoXP)){
    if(NoXP < 3) stop("NoXP is too small")
    NoXP <- as.integer(NoXP)
  }
  n <- length(x)
  if(setDelta <- is.null(delta)){
    if(n <= 100) delta <- 0
    else delta <- diff(range(x))/100  ## if the values in x are uniformly
    ## scattered over the range, the local weighted regression would be
    ## carried out at approximately 100 points.
  }
  ## weights
  if(is.null(weight)) weight <- rep(1,n)
  ## initial fit
  if(is.null(init)){
    rw <- rep(1, n)
    MAXIT <- c(maxit[1],sum(maxit)) + 1
  }
  else{
    if(length(init) != n)
      stop("Must have same number of initial values as observations")
    Resid <- (y - init)*sqrt(weight)
    scale <- Scale(Resid)
    if(scale == 0)
      status <- stop("could not compute scale of residuals for iter = 0")
    else{
      if(maxit[1]==0){
	u <- abs(Resid/(scale*b))
	rw <- ifelse(Resid < 0, 1, ifelse(u > 1, 0, ((1 + u) * (1 - u))^2))
      }
      else{
	u <- abs(Resid/(scale*b))
	rw <- ifelse(u > 1, 0, ((1 + u) * (1 - u))^2)
      }
    }
    MAXIT <- c(maxit[1],sum(maxit))
  }
  if(SORT){
    h <- order(x)
    x <- x[h]
    y <- y[h]
    weight <- weight[h]
    rw <- rw[h]
  }
  ## load C function if necessary (Complile by Rcmd SHLIB lwreg.c)
  if(!is.loaded("lwreg")){
    dyn.load("lwreg.dll")
    on.exit(dyn.unload("lwreg.dll"))
  }
  fit <- NULL
##
## "robust iteration"  (irlls)
  for(iiter in 1:MAXIT[2]) {
    if(DOT){
      ok <- rw > 10^(-6)
      xx <- x[ok]
      wweight <- weight[ok]*rw[ok]
      yy <- y[ok]
      nn <- sum(ok)
      if(setDelta){
        if(nn <= 100) delta <- 0
        else delta <- diff(range(x))/100  ## cf at top of this procedure
      }
      qub <- if(is.null(NoXP)) max(floor(span*nn), 2)
             else NoXP
      fit <- .C('lwreg',
                PACKAGE="IDPmisc",
		xx=as.double(xx),
		as.double(yy),
		as.integer(nn),
		as.integer(qub),
		as.double(delta),
		as.double(wweight),
		yfit=double(nn))
    }
    else{
      qub <- if(is.null(NoXP)) max(floor(span*n), 2)
             else NoXP
      wweight <- weight*rw
      fit <- .C('lwreg',
                PACKAGE="IDPmisc",
		xx=as.double(x),
		as.double(y),
		as.integer(n),
		as.integer(qub),
		as.double(delta),
		as.double(wweight),
		yfit=double(n))
    }
    ## calculation of weights for irls
    if(iiter < MAXIT[2]){
      ## compute robustness weights except last iteration
      Resid <- (y - approx(x=fit$xx, y=fit$yfit, xout=x, rule=2)$y)*sqrt(weight)
      scale <- Scale(Resid)
      if(scale == 0)
	status <- stop(paste("could not compute scale of residuals for iter =",
			     iiter))
      else{
 	if(iiter < MAXIT[1]){
 	  u <- abs(Resid/(scale*b))
 	  rw <- ifelse(Resid < 0, 1, ifelse(u > 1, 0, ((1 + u) * (1 - u))^2))
 	}
 	if(iiter >= MAXIT[1]){
 	  u <- abs(Resid/(scale*b))
 	  rw <- ifelse(u > 1, 0, ((1 + u) * (1 - u))^2)
	}
      }
    }  ## end of irlls
  }
##
  list(x=x, y=y, fit=approx(x=fit$xx, y=fit$yfit, xout=x, rule=2)$y,
       rw=rw, sigma=scale, span=span, NoXP=NoXP, maxit=maxit, b=b,
       weight=weight, Scale=Scale, delta=delta, SORT=SORT, DOT=DOT)
} ## rfbaseline

