"ToeplitzBlock" <- 
  function(res,lag.max,Kernel=FALSE)
  {
    res <- as.matrix(res)
    k <- NCOL(res)
    n <- NROW(res)
    m <- lag.max+1
    out <- matrix(numeric(k*m*k*(lag.max+1)),nrow=k*m,ncol=k*m)
    if (Kernel==FALSE){
      Accmat <- stats::acf(res, lag.max = lag.max, plot = FALSE, type = "covariance")$acf
      inveseC0 <- solve(Accmat[1,,])
      L <- t(chol(inveseC0))
      for (i in 0:lag.max)
        for (j in i:lag.max){
          out[(j*k+1):(k*(j+1)),(i*k+1):(k*(i+1))] <- crossprod(t(crossprod(L,t(Accmat[j-i+1,,]))),L)
          out[(i*k+1):(k*(i+1)),(j*k+1):(k*(j+1))] <- crossprod(t(crossprod(L,Accmat[j-i+1,,])),L)
        }
    }
    else 
    {
      accr<-stats::acf(res,lag.max,plot=FALSE,type="correlation")$acf

      x<-(0:lag.max)/lag.max
      weights<-ifelse( x > 1, 0, ifelse( x < 0.5, 1 - 6 *  x^2 + 6 * abs( x)^3, 2 * (1 - abs( x))^3))
      weightaccr<-weights*accr
      out<-toeplitz(as.vector(weightaccr))
    }
    return(out)
  }
