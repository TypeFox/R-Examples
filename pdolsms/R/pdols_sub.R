# @title Procedures for Long-run Estimators
# @description 
# se =  plong(x,p);  Paramteric method for long-run estimator with fixed lag length.                               @
# se =  pdlong(x,p): Parametric method for long-run estimator with automatic lag selection based on Hall           @
# se =  nwest(x,p):  Non-parametric method for long-run estimator Newey-West Batlett Kernel estimator             @
# @param x (Tx1) vector                                                   
# @param p lag length                                                    
# @return se long-run standard error                                               @


plong <- function(x, p) {
  if (NCOL(x) > 1) break("x must be a Tx1 vector")
  t <- NROW(x)
  x1 <- matrix(0, nrow = t - p, ncol = p + 1)
  j <- 0
  while (j <= p) {
    x1[, j+1] <- x[(p-j+1):(NROW(x) - j), ]
    j <- j + 1
  }
  yx <- x1[, 1, drop = FALSE]
  xx <- x1[, 2:NCOL(x1), drop = FALSE]
  b <- solve(crossprod(xx)) %*% crossprod(xx, yx)
  resi <- yx - xx %*% b
  s_s <- crossprod(resi)/(NROW(x) - NROW(b))
  lambda <- s_s / (1 - apply(b, 2, sum)) ^ 2
  
  return(as.numeric(lambda))
  
}


pdlong <- function(x, p) {
  if (NCOL(x) > 1) break("x must be a Tx1 vector")
  t <- NROW(x)
  for (p in p:0) {
    ip <- p
    if(p == 0) {
      lambda <- sqrt(crossprod(x)/(NROW(x) - 1))
      return(as.numeric(lambda))    
    }
    
    x1 <- matrix(0, nrow = t - ip, ncol = ip + 1)
    for (j in 0:ip) {
      x1[, j+1] <- x[(ip-j+1):(NROW(x) - j), , drop = FALSE]
    }  
    
    y <- x1[, 1]
    xx <- x1[, 2:NCOL(x1)]
    b <- solve(crossprod(xx)) %*% crossprod(xx, y)
    resi <- y - xx %*% b
    s_s <- as.numeric(crossprod(resi)/(NROW(y) - NROW(b)))
    ss <- diag(s_s * solve(crossprod(xx)))
    tl <- abs(b[NROW(b)]) / sqrt(ss[NROW(ss)])
    
    if (tl > 1.96) {
      lambda <- s_s/(1-apply(b, 2, sum))^2
      return(as.numeric(lambda))
    }
    
  }
}



nwest <- function(x, inlag) {
  if (NCOL(x) > 1) break("x must be a Tx1 vector")
  n <- NROW(x)
  xom <- crossprod(x)/n
  cat(xom, "\n")
  xomm <- 0
  
  for (i1 in 1:inlag) {
    temp11 <- apply((x[(1+i1):(n-1), , drop = FALSE] * 
                       x[1:(n-1-i1), , drop = FALSE]), 
                    2, sum)
    
    xomm <- (1-(i1/(inlag+1)))*(temp11/n) + xomm
  }
  
  xom <- xom + 2 * xomm
  
  return(as.numeric(xom))
}


sauto <- function(x, p) {
  if (NCOL(x) > 1) break("x must be a Tx1 matrix")
  z <- matrix(0, NROW(x) - p, p + 1)
  
  for(i in 0:p) {
    z[, i + 1] <- x[(p - i + 1):(NROW(x) - i), ]
  }
  
  return(z) 
  # First column = y1 and rest of columns = x1
  
}


qspw <- function(x) {
  if (NCOL(x) > 1) break("x must be a Tx1 vector")
  t <- NROW(x)
  z <- sauto(x, 1)
  uy <- z[, 1]
  ux <- z[, 2:NCOL(z)]
  
  b1 <- solve(crossprod(ux)) %*% crossprod(ux, uy)
  if (b1 >= 1) {
    b1 <- 0.999
  }
  
  r1 <- uy - ux * b1
  a2 <- (4 * b1^2)/ ((1-b1)^4)
  # browser()
  mq <- floor(1.3221 * (a2 * t) ^ (1/5))
  
  if (mq == 0) {
    mq <- 1
  }
  if (mq >= (t-2)) {
    mq <- t - 4
  }
  
  ome <- matrix(0, nrow = t-2, ncol = 1)
  ome0 <- crossprod(r1)/(t - 1)
  kq <- ome
  for (i in 1:(t-1-2)) {
    ome[i, 1] <- crossprod(r1[(i + 1):(t - 1)], r1[1:(t-1-i)])/(t-1)
    kq[i, 1] <- (25/ (12 * pi^2 * (i/mq)^2)) *
      ((sin( 6*pi*(i/mq)/5 )/(6*pi*(i/mq)/5) )
        - cos(6*pi*(i/mq)/5))
  }
  
  ssq <- ome0 + 2 * apply(ome * kq, 2, sum)

  ssq <- ssq / ((1-b1)^2)
  # cat(ssq, "\n")
  return(as.numeric(ssq))
  
}


