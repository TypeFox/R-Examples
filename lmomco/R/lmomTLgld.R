"lmomTLgld" <-
function(para, nmom=6, trim=1, leftrim=NULL, rightrim=NULL, tau34=FALSE) {
  L <- R <- vector(mode="numeric", length=nmom)

  if(nmom < 2) {
    warning("Number of L-moments is less than 2")
    return()
  }
  if(! is.null(trim) && trim < 0) {
    warning("Trimming value is less than 0")
    return()
  }
  if(! is.null(leftrim) && leftrim < 0) {
    warning("Left rimming value is less than 0")
    return()
  }
  if(! is.null(rightrim) && rightrim < 0) {
    warning("Right trimming value is less than 0")
    return()
  }
  if(is.null(trim) && is.null(leftrim) && is.null(rightrim)) {
    trim <- 0
  }
  if(is.null(leftrim))  leftrim  <- trim
  if(is.null(rightrim)) rightrim <- trim

  t1 <- leftrim
  t2 <- rightrim

  E <- para$para[1]
  A <- para$para[2]
  K <- para$para[3]
  H <- para$para[4]

  if(K <= -(1+leftrim)) {
     warning("Parameter k is too small for the leftrim level: k > -(1+leftrim), still yet to check the other three")
     return()
  }
  if(H <= -(1+rightrim)) {
     warning("Parameter h is too small for the rightrim level: k > -(1+rightrim), still yet to check the other three")
     return()
  }


  #if((t1 + t2) == 0) {
     if(! are.pargld.valid(para)) return()
     attributes(para$para) <- NULL
  #}


  if(tau34) {
    for(r in 3:4) {
      the.sum <- 0
      for(j in 0:(r-1)) {
        sig  <- (-1)^j
        tmpA <- choose(r-1, j)*choose(r+t1+t2-1,r+t1-j-1)
        a <- lgamma(K+r+t1-j)
        b <- lgamma(t2+j+1)
        c <- lgamma(K+r+t1+t2+1)
        d <- lgamma(r+t1-j)
        e <- lgamma(H+t2+j+1)
        f <- lgamma(H+r+t1+t2+1)
        tmpB <- exp(a+b-c) - exp(d+e-f)
        the.sum <- the.sum + sig*tmpA*tmpB
      }
      L[r] <- A * (r+t1+t2) * the.sum / r
    }
    R[1] <- NA
    R[2] <- NA
    for(r in 3:4) R[r] <- L[r]/L[2]
    return(R[3], R[4])
  } else {
  for(r in 1:nmom) {
    if(r > 1) E <- 0
    the.sum <- 0
    for(j in 0:(r-1)) {
      sig  <- (-1)^j
      tmpA <- choose(r-1, j)*choose(r+t1+t2-1,r+t1-j-1)
      a <- lgamma(K+r+t1-j)
      b <- lgamma(t2+j+1)
      c <- lgamma(K+r+t1+t2+1)
      d <- lgamma(r+t1-j)
      e <- lgamma(H+t2+j+1)
      f <- lgamma(H+r+t1+t2+1)
      tmpB <- exp(a+b-c) - exp(d+e-f)
      the.sum <- the.sum + sig*tmpA*tmpB
    }
    L[r] <- E + A * (r+t1+t2) * the.sum/ r
  }
  R[1] <- NA
  R[2] <- L[2]/L[1]
  if(nmom >= 3) for(r in 3:nmom) R[r] <- L[r]/L[2]
  if(t1 == t2) {
     trim <- t1
  } else {
     trim <- NULL
  }
  z <- list(lambdas = L, ratios = R,
            trim=trim, leftrim=leftrim, rightrim=rightrim,
            source = "lmomTLgld")
  return(z)
  }
}

