harvtest <- function(formula, order.by=NULL, data=list())
{
  dname <- paste(deparse(substitute(formula)))

  if(!inherits(formula, "formula")) {
    X <- if(is.matrix(formula$x))
           formula$x
         else model.matrix(terms(formula), model.frame(formula))
    y <- if(is.vector(formula$y))
           formula$y
         else model.response(model.frame(formula))
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
  }  

  k <- ncol(X)
  n <- nrow(X)

  rec.res <- function(X, y)
  {
      n <- nrow(X)
      q <- ncol(X)
      w <- rep(0,(n-q))
      Xr1 <- X[1:q,,drop = FALSE]
      xr <- as.vector(X[q+1,])
      X1 <- chol2inv(qr.R(qr(Xr1)))
      fr <- as.vector((1 + (t(xr) %*% X1 %*% xr)))
      betar <- X1 %*%t(Xr1)%*% y[1:q]
      w[1] <- (y[q+1] - t(xr) %*% betar)/sqrt(fr)

      for(r in ((q+2):n))
      {
          X1 <- X1 - (X1 %*% outer(xr, xr) %*% X1)/fr
  	  betar <- betar + X1 %*% xr * w[r-q-1]*sqrt(fr)
  	  xr <- as.vector(X[r,])
          fr <- as.vector((1 + (t(xr) %*% X1 %*% xr)))
          w[r-q] <- (y[r] - t(xr) %*% betar)/sqrt(fr)
      }
      return(w)
  }

  if(!is.null(order.by))
  {
    if(inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[,ncol(z)])
    } else {
      z <- order.by
    }
    X <- as.matrix(X[order(z),])
    y <- y[order(z)]
  }

  resr <- rec.res(X,y)
  sigma <- sqrt(var(resr)*(length(resr)-1)/(n-k-1))
  resr <- resr / sigma
  harv <- abs(sum(resr)/sqrt(n-k))/sqrt(var(resr))
  names(harv) <- "HC"
  df <- n-k-1
  names(df) <- "df"
  RVAL <- list(statistic = harv,
      parameter = df,
      method = "Harvey-Collier test",
      p.value= 2 * pt(harv, n-k-1,lower.tail=FALSE),
      data.name=dname)

  class(RVAL) <- "htest"
  return(RVAL)
}

