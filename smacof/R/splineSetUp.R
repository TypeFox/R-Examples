splineSetUp <- function(x,allKnots,degree){
  # Sets up a monotone spline base.
  #   AllKnots: the number of interior and exterior knots
  #   degree:    degree of the spline
  knotSeq <- createKnots(x,allKnots,degree)
  base    <- JSpline(x,knotSeq,degree,allKnots)
  # Remove knot duplications at extrema
  knotSeq <- knotSeq[(degree-1)+(1:allKnots)]
  return(list(base = base, knotSeq = knotSeq))
}


createKnots <- function(x,AllKnots,degree){
  # Creates the knot sequence.
  #   AllKnots: the number of interior and exterior knots
  #   degree:    degree of the spline
  #   
  #
  n    <- length(x)
  Temp <- sort(x)
  m0   <- degree
  m1   <- m0 + 1
  m2   <- degree + AllKnots - 2
  m3   <- m2 + 1
  m4   <- m3 + degree - 1
  tt   <- rep(0,m4)
  tt[1:m0]<- Temp[1]
  tt[(m3-1)+(1:(m4-m3+1))] <- Temp[n]
  xdegree <- m3 - m0
  prop   <- ((m1 - m0 - 1) + (1:(m2-m1+1)))/xdegree
  xpos  <- prop*(n+1)
  npos  <- floor(xpos)
  xval  <- Temp[npos]+(xpos-npos)*(Temp[npos+1]-Temp[npos])
  tt[(m1-1)+(1:(m2-m1+1))] <- xval
  knotSeq <- tt    
  return(knotSeq)
}

ImSpline <- function(knotSeq,datum,degree,left){ 
  isw   <- 0
  ncoef <- length(knotSeq) - degree
  wk    <- rep(1,ncoef)
  dr    <- wk
  dl    <- wk
  spli  <- wk
  j     <- 1
  while (isw==0){
    isw <- degree==1
    if (j>=degree) isw=1
    jpl <- j+1
    dr[j] <- knotSeq[left+j] - datum
    dl[j] <- datum - knotSeq[left-j+1]
    s <- 0
    for (i in 1:j){
      z     <- wk[i]/(dr[i]+dl[jpl-i])
      wk[i] <- s+dr[i]*z
      s     <- dl[jpl-i]*z
    }
    wk[jpl] <- s
    if (isw==1) break 
    j <- jpl
  }
  s <- 0
  for (i in 1:degree){
    n <- degree - i + 1
    s <- s + wk[n+1]
    spli[n] <- s
  }
  return(spli)
}

JSpline <-  function(x,knotSeq,degree,AllKnots){
  # Compute I-spline basis (base) for variable x
  n     <- length(x)
  ncoef <- degree + AllKnots - 2
  base  <- matrix(0,n,ncoef)
  for (i in 1:n){
    datum <- x[i]
    left  <- degree
    if (left <= ncoef - 1){
      while(left <= ncoef - 1){
        if (knotSeq[left+1] > datum) break
        left <- left + 1
      }
    }
    #sum((degree:(ncoef-1))<datum);
    spli <- ImSpline(knotSeq, datum, degree, left) 
    if (left > degree) {
      ind  <- 1:(left-degree)
      base[i,ind] <- base[i,ind] + 1
    }
    ind  <- (left-degree+1):left 
    base[i,ind] <- base[i,ind] + spli[1:length(ind)]
  }
  return(base)
}
