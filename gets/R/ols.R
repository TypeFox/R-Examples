ols <-
function(y, x, tol = 1e-07 , LAPACK = FALSE,
  method=1, user.fun=NULL, user.options=NULL)
{
  ##user-specified:
  if(method==0){
    user.options <- c(list(y=y, x=x), user.options)
    out <- do.call(user.fun, user.options)
  }

  ##fastest (usually only for estimates):
  if(method==1){
    out <- list()
    qx <- qr(x, tol, LAPACK = LAPACK)
    out <- c(out, qx)
    out$coefficients <- as.vector(solve.qr(qx, y, tol = tol))
  }

  ##second fastest:
  if(method==2){
    tmp <- crossprod(x)
    out <- qr(tmp, tol, LAPACK=LAPACK)
    out$xtxinv <- solve.qr(out, tol=tol, LAPACK=LAPACK)
    out$xtx <- tmp
    out$xty <- crossprod(x,y)
    out$coefficients <- as.vector(out$xtxinv%*%out$xty)
    out$fitted <- as.vector(x %*% cbind(out$coefficients))
    out$residuals <- y - out$fitted
  }

  ##ordinary vcov:
  if(method==3){
    tmp <- crossprod(x)
    out <- qr(tmp, tol, LAPACK=LAPACK)
    out$xtxinv <- solve.qr(out, tol=tol, LAPACK=LAPACK)
    out$xtx <- tmp
    out$xty <- crossprod(x,y)
    out$coefficients <- as.vector(out$xtxinv%*%out$xty)
    out$fitted <- as.vector(x %*% cbind(out$coefficients))
    out$residuals <- y - out$fitted
    out$resids2 <- out$residuals^2
    out$rss <- sum(out$resids2)
    out$df <- length(y) - NCOL(x)
    out$sigma2 <- out$rss/out$df
    out$vcov <- out$sigma2*out$xtxinv
  }

  ##white vcov:
  if(method==4){
    tmp <- crossprod(x)
    out <- qr(tmp, tol, LAPACK=LAPACK)
    out$xtxinv <- solve.qr(out, tol=tol, LAPACK=LAPACK)
    out$xtx <- tmp
    out$xty <- crossprod(x,y)
    out$coefficients <- as.vector(out$xtxinv%*%out$xty)
    out$fitted <- as.vector(x %*% cbind(out$coefficients))
    out$residuals <- y - out$fitted
    out$resids2 <- out$residuals^2
    out$rss <- sum(out$resids2)
    out$df <- length(y) - NCOL(x)
    out$sigma2 <- out$rss/out$df
    out$omegahat <- crossprod(x, x*out$resids2)
    out$vcov <- out$xtxinv %*% out$omegahat %*% out$xtxinv
  }

  ##newey-west vcov:
  if(method==5){
    tmp <- crossprod(x)
    out <- qr(tmp, tol, LAPACK=LAPACK)
    out$xtxinv <- solve.qr(out, tol=tol, LAPACK=LAPACK)
    out$xtx <- tmp
    out$xty <- crossprod(x,y)
    out$coefficients <- as.vector(out$xtxinv%*%out$xty)
    out$fitted <- as.vector(x %*% cbind(out$coefficients))
    out$residuals <- y - out$fitted
    out$resids2 <- out$residuals^2
    out$rss <- sum(out$resids2)
    out$df <- length(y) - NCOL(x)
    out$sigma2 <- out$rss/out$df

    y.n <- length(y)
    iL <- round(y.n^(1/4), digits=0)
    vW <- 1 - 1:iL/(iL+1)
    vWsqrt <- sqrt(vW)
    mXadj <- out$residuals*x
    mS0 <- crossprod(mXadj)

    mSum <- 0
    for(l in 1:iL){
      mXadjw <- mXadj*vWsqrt[l]
      mXadjwNo1 <- mXadjw[-c(1:l),]
      mXadjwNo2 <- mXadjw[-c(c(y.n-l+1):y.n),]
      mSum <- mSum + crossprod(mXadjwNo1, mXadjwNo2) + crossprod(mXadjwNo2, mXadjwNo1)
    }

    out$omegahat <- mS0 + mSum
    out$vcov <- out$xtxinv %*% out$omegahat %*% out$xtxinv
  }

  ##result:
  return(out)

}
