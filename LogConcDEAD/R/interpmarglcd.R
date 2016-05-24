##This function is internal and computes the integral of
## exp(b^T x - beta) over
## T \cup L (w.r.t. Lebesgue measure on an appropriate subspace)
## where T is the simplex with vertices given by verts
## and L is the subspace X[,keep] = lambda

'integratemarg' <- function(keep, lambda, b, beta, verts, d=2, eps=10^-6) {
  ## Find vertices are to the left, right or on the current vertex
  vals <- verts[,keep] - lambda
  left <- vals < -eps
  right <- vals > eps
  on <- abs(vals) <= eps

  ans <- 0
  new <- NULL
  new <- rbind(new,verts[on,])
  for (i in which(left)) {
    for (j in which(right)) {
      new <- rbind(new, verts[i,] + (lambda - verts[i,keep])/(verts[i,keep] - verts[j,keep]) *(verts[i,] - verts[j,]))
    }
  }
  if(nrow(new)==d) {
    tmp <- new[,-keep,drop=FALSE]
    y <- new%*%b - beta
    ans <- ans + JAD(y)%*%abs(det(tmp[-1,,drop=FALSE]-rep(1,d-1)%*%tmp[1,,drop=FALSE]))
  }
  else if(nrow(new)>d) {
    tmp <- new[,-keep,drop=FALSE]
    trig <- matrix(delaunaynew(tmp),ncol=d)
    for (i in 1:nrow(trig)) {
      y <- new[trig[i,],]%*%b - beta
      ans <- ans + JAD(y)%*%abs(det(tmp[trig[i,-1],,drop=FALSE]-rep(1,d-1)%*%tmp[trig[i,1],,drop=FALSE]))
    }
  }
  return(ans)
}


##Computes one-dimensional marginals!

'interpmarglcd' <- function(lcd, marg=1, gridlen=100) {
  if(class(lcd) != "LogConcDEAD") {
    stop("error: lcd must be of class LogConcDEAD")
  }
  x <- lcd$x
  d <- ncol(x)
  if(!is.element(marg,1:d)) {
    stop("error: marg must be one of 1 ... ",d)
  }
  triang <- lcd$triang
  xo <- seq(min(x[,marg]),max(x[,marg]),len=gridlen)
  val <- rep(0,gridlen)
  for (j in 1:nrow(triang)) {
    for (i in which(xo > min(x[triang[j,],marg]) & xo < max(x[triang[j,],marg]))) {
      val[i] <- val[i] + integratemarg(marg, xo[i], lcd$b[j,], lcd$beta[j], lcd$x[lcd$triang[j,],],d=d)
    }
  }
  return(list(xo=xo,marg=val))
}


'dmarglcd' <- function(x=0, lcd,  marg=1) {
  if(class(lcd) != "LogConcDEAD") {
    stop("error: lcd must be of class LogConcDEAD")
  }
  d <- ncol(lcd$x)
  if(!is.element(marg, 1:d)) {
    stop("error: marg must be one of 1 ... ",d)
  }
  if(!is.numeric(x)) {
    stop("x must be numeric")
  }
  triang <- lcd$triang
  ansvec <- NULL
  for (i in 1:length(x)) {
    ans <- 0
    for (j in 1:nrow(triang)) {
      ans <- ans + integratemarg(marg, x[i], lcd$b[j,], lcd$beta[j], lcd$x[lcd$triang[j,],],d=d)
    }
    ansvec <- c(ansvec,ans)
  }
  return(ansvec)
}
