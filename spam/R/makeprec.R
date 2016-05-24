# This is file ../spam/R/makeprec.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     



"precmat.GMRFreglat" <- function(n,m, par=.1, model = "m1p1",  eps = .Spam$eps){
  if((n<2)|(m<2))
    stop("n and m need to be >1")
  dims <- c(n,m)
  if (model=="m1p1"){
    x <- numeric(dims[1])
    x[1:2] <- c(1,-par[1])
    y <- numeric(prod(dims))
    y[dims[1]+1] <- -par[1] 
    return( kronecker(diag.spam(dims[2]), toeplitz.spam(x,eps=eps))+toeplitz.spam(y,eps=eps))
  }
  
  if (model=="m1p2"){
    x <- numeric(dims[1])
    x[1:2] <- c(1,-par[1])
    y <- numeric(prod(dims))
    y[dims[1]+1] <- -par[2] 
    return( kronecker(diag.spam(dims[2]), toeplitz.spam(x,eps=eps))+toeplitz.spam(y,eps=eps))
  }
  
  if (model=="m2p3"){
    x <- numeric(dims[1])
    x[1:2] <- c(1,-par[1])
    
    y <- numeric(dims[1])
    y[1:2] <- c(-par[2],-par[3])

    z <- numeric(dims[2])
    z[2] <- 1


    p1 <- kronecker( diag.spam(dims[2]), toeplitz.spam(x,eps=eps))
    p2 <- kronecker( toeplitz.spam(z,eps=eps), toeplitz.spam(y,eps=eps))
    return( p1 + p2)
  }
  if (model=="m2p4"){
    x <- numeric(dims[1])
    x[1:2] <- c(1,-par[1])
    
    y <- numeric(dims[1])
    y[1:2] <- c(-par[2],-par[3])
    w <- numeric(dims[1])
    w[1:2] <- c(-par[2],-par[4])

    z <- numeric(dims[2])
    z[2] <- 1


    p1 <- kronecker( diag.spam(dims[2]), toeplitz.spam(x,eps=eps))
    p2 <- kronecker( toeplitz.spam(z,rep(0,dims[2]),eps=eps), toeplitz.spam(y,w,eps=eps))
    p3 <- kronecker( toeplitz.spam(rep(0,dims[2]),z,eps=eps), toeplitz.spam(w,y,eps=eps))
    return( p1 + p2 + p3)
  }
  stop("Model not implemented yet!")

}
