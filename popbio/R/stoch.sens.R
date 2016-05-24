stoch.sens <- function(A, tlimit = 100){

# june 30, 2015 added option for list of matrices..

if(!is.list(A)  && !is.matrix(A[[1]])) stop("A should be a list of matrices")
k <- ncol(A[[1]])

A <- sample(A, tlimit, replace=TRUE)


# A translation of Caswell's (2001) Matlab code fragment
  tlimit <- length(A)
  wvec <- rep(1/k,k)
  w <- cbind(wvec)

  # generate sequence of structure vectors

  r <- rep(0,tlimit)

  for(i in 1:tlimit){
    a <- A[[i]]
    wvec <- a %*% wvec
    r[i] <- sum(wvec)
    wvec <- wvec/r[i]
    w <- cbind(w,wvec)
  }

  # specifiy initial reproductive value vector

  vvec <- rep(1/k,k)
  v <- cbind(vvec)

  for(i in rev(1:tlimit)){

    a <- A[[i]]
    vvec <- vvec%*%a
    v <- cbind(t(vvec),v)
  }

  sensmat <- matrix(0,nrow=k,ncol=k)
  elasmat <- matrix(0,nrow=k,ncol=k)

  for(i in 1:tlimit){
    # for some reason, need the as.numeric() to get the division by
    # scalar to work 
    sensmat <- sensmat+((v[,i+1]%*%t(w[,i])) /
                        as.numeric(r[i]*t(v[,i+1])%*%w[,i+1]))
    a <- A[[i]]
    elasmat <- elasmat+((v[,i+1]%*%t(w[,i])*a) /
                        as.numeric((r[i]*t(v[,i+1])%*%w[,i+1])))
  }

  sensmat <- sensmat/tlimit
  elasmat <- elasmat/tlimit

  out <- list(sensitivities=sensmat, elasticities=elasmat)

  out
}
