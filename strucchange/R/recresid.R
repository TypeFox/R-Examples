recresid <- function(x, ...)
{
  UseMethod("recresid")
}

recresid.formula <- function(formula, data = list(), ...)
{
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    modelterms <- terms(formula, data = data)
    X <- model.matrix(modelterms, data = data)
    rr <- recresid(X, y, ...)
    return(rr)
}

recresid.lm <- function(x, data = list(), ...)
{
    X <- if(is.matrix(x$x)) x$x else model.matrix(terms(x), model.frame(x))
    y <- if(is.vector(x$y)) x$y else model.response(model.frame(x))
    rr <- recresid(X, y, ...)
    return(rr)
}

recresid.default <- function(x, y, start = ncol(x) + 1, end = nrow(x),
  tol = sqrt(.Machine$double.eps)/ncol(x), ...)
{
  ## checks and data dimensions
  stopifnot(start > ncol(x) & start <= nrow(x))
  stopifnot(end >= start & end <= nrow(x))
  n <- end
  q <- start - 1
  k <- ncol(x)
  rval <- rep(0, n - q)

  ## convenience function to replace NAs with 0s in coefs
  coef0 <- function(obj) {
    cf <- obj$coefficients
    ifelse(is.na(cf), 0, cf)
  }
  Xinv0 <- function(obj) {
    qr <- obj$qr
    rval <- matrix(0, ncol = k, nrow = k)
    wi <- qr$pivot[1:qr$rank]
    rval[wi,wi] <- chol2inv(qr$qr[1:qr$rank, 1:qr$rank, drop = FALSE])
    rval
  }

  ## initialize recursion
  y1 <- y[1:q]
  fm <- lm.fit(x[1:q, , drop = FALSE], y1)
  X1 <- Xinv0(fm)	 
  betar <- coef0(fm)
  xr <- as.vector(x[q+1,])
  fr <- as.vector((1 + (t(xr) %*% X1 %*% xr)))
  rval[1] <- (y[q+1] - t(xr) %*% betar)/sqrt(fr)

  ## check recursion agains full QR decomposition?
  check <- TRUE

  if((q+1) < n)
  {
    for(r in ((q+2):n))
    {
        ## check for NAs in coefficients
	nona <- all(!is.na(fm$coefficients))
    
	## recursion formula
        X1 <- X1 - (X1 %*% outer(xr, xr) %*% X1)/fr
        betar <- betar + X1 %*% xr * rval[r-q-1] * sqrt(fr)

	## full QR decomposition
	if(check) {
	  y1 <- y[1:(r-1)]
	  fm <- lm.fit(x[1:(r-1), , drop = FALSE], y1)
	  nona <- nona & all(!is.na(betar)) & all(!is.na(fm$coefficients))
	  ## keep checking?
	  if(nona && isTRUE(all.equal(as.vector(fm$coefficients), as.vector(betar), tol = tol))) check <- FALSE 
	  X1 <- Xinv0(fm)
          betar <- coef0(fm)
        }
        
        ## residual
        xr <- as.vector(x[r,])
	fr <- as.vector((1 + (t(xr) %*% X1 %*% xr)))
	rval[r-q] <- (y[r] - sum(xr * betar, na.rm = TRUE))/sqrt(fr)
    }
  }
  return(rval)
}

