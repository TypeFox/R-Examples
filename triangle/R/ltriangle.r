################################################################################
#
# Program:   ltriangle.R
# Purpose:   To draw from the log10 triangle distribution
# Author:    Rob Carnell
# Date:      June 06
#
# Variables
#   used the same naming conventions as other R distributions (r,q,p,d)
#   n = number of values to return
#   a = left triangle endpoint,
#   b = right triangle endpoint
#   c = distribution mode
#
################################################################################

rltriangle <- function(n=1, a=1, b=100, c=10^((log10(a)+log10(b))/2), logbase=10)
{
  stopifnot(length(n) == 1)
  if(n < 1 | is.na(n)) stop(paste("invalid argument: n =", n))
  n <- floor(n)
  if(any(is.na(c(a,b,c)))) return(rep(NaN, times=n)) # to match behavior of runif
  if(any(a > c | b < c)) return(rep(NaN, times=n)) # to match behavior of runif
  if(any(is.infinite(c(a, b, c)))) return(rep(NaN, times=n))
  if(any(c(a,b,c) == 0)) return(rep(-Inf, times=n))
  if(any(c(a,b,c) < 0)) return(rep(NaN, times=n))
  
  lp <- runif(n)

  stopifnot(length(logbase) == 1)
  if (logbase==10)
  {
    la <- log10(a)
    lb <- log10(b)
    lc <- log10(c)
  } else
  {
    la <- log(a)/log(logbase)
    lb <- log(b)/log(logbase)
    lc <- log(c)/log(logbase)
  }
  
  if(a != c)
  {
    # if a = c then i is always true
    i <- which((la + sqrt(lp * (lb - la)*(lc - la))) <= lc)
    j <- which((lb - sqrt((1 - lp) * (lb - la) * (lb - lc))) > lc)
  } else
  {
    i <- which((la + sqrt(lp * (lb - la)*(lc - la))) < lc)
    j <- which((lb - sqrt((1 - lp) * (lb - la) * (lb - lc))) >= lc)
  }

  if(length(i) != 0)
    lp[i] <- la + sqrt(lp[i] * (lb - la) * (lc - la))
  if(length(j) != 0)
    lp[j] <- lb - sqrt((1 - lp[j]) * (lb - la) * (lb - lc))
    
  p <- logbase^lp

  return(p)
}

dltriangle <- function(x, a=1, b=100, c=10^((log10(a)+log10(b))/2), logbase=10) {
  x1 <- length(x)
  a1 <- length(a)
  b1 <- length(b)
  c1 <- length(c)

  stopifnot(length(logbase) == 1)
  if (logbase==10)
  {
    la <- log10(a)
    lb <- log10(b)
    lc <- log10(c)
    lx <- log10(x)
  } else
  {
    la <- log(a)/log(logbase)
    lb <- log(b)/log(logbase)
    lc <- log(c)/log(logbase)
    lx <- log(x)/log(logbase)
  }

  dTest <- function(X){
    if(any(is.na(X)))
    {
      # is.na is TRUE for NA, NaN, and FALSE
      if(any(is.nan(X))) return(NaN) # to conform to qunif
      else return(NA) # to conform to qunif
    } else if(X[2] > X[4] | X[3] < X[4] | (X[1]==X[2] & X[2]==X[4]))
    {
      warning("values required to be  a <= c <= b (at least one strict inequality)")
      return(NaN) # to conform to behavior of qunif
    } else if(any(is.infinite(X[2:4])))
    {
      return(NaN)
    } else if(X[1] <= X[2])
    {
      return(0)
    } else if(X[2] != X[4] & X[1] < X[4])
    {
      return(2*(X[1] - X[2]) / (X[3] - X[2]) / (X[4] - X[2]))
    } else if(X[4] != X[3] & X[1] >= X[4] & X[1] < X[3])
    {
      return(2*(X[3] - X[1]) / (X[3] - X[2]) / (X[3] - X[4]))
    } else if(X[1] >= X[3]){
      return(0)
    }
  }

  k <- max(x1, a1, b1, c1)
  if(k==1) return(dTest(c(lx, la, lb, lc)))

  params <- matrix(nrow=k, ncol=4)
  tryCatch(
  {
    params[,1] <- lx
    params[,2] <- la
    params[,3] <- lb
    params[,4] <- lc
  }, error = function(X) {
    stop(paste(" -- Argument Lengths: length of x = ", x1,
                ", a = ", a1, ", b = ", b1, ", c = ", c1, " -- ", X, sep=""))
  })

  return(apply(params, 1, dTest))
}

pltriangle <- function(q, a=1, b=100, c=10^((log10(a)+log10(b))/2), logbase=10)
{
  q1 <- length(q)
  a1 <- length(a)
  b1 <- length(b)
  c1 <- length(c)

  stopifnot(length(logbase) == 1)
  if (logbase==10)
  {
    la <- log10(a)
    lb <- log10(b)
    lc <- log10(c)
    lq <- log10(q)
  } else
  {
    la <- log(a)/log(logbase)
    lb <- log(b)/log(logbase)
    lc <- log(c)/log(logbase)
    lq <- log(q)/log(logbase)
  }

  pTest <- function(X)
  {
    if(any(is.na(X))){ # is.na is TRUE for NA, NaN, and FALSE
      if(any(is.nan(X))) return(NaN) # to conform to qunif
      else return(NA) # to conform to qunif
    } else if(X[2] > X[4] | X[3] < X[4] | (X[1]==X[2] & X[2]==X[4]))
    {
      warning("values required to be  a <= c <= b (at least one strict inequality)")
      return(NaN) # to conform to behavior of qunif
    } else if(any(is.infinite(X[2:4])))
    {
      return(NaN)
    } else if(X[1] <= X[2])
    {
      return(0)
    } else if(X[2] != X[4] & X[1] < X[4])
    {
      return((X[1] - X[2])^2 / (X[3] - X[2]) / (X[4] - X[2]))
    } else if(X[4] != X[3] & X[1] >= X[4] & X[1] < X[3])
    {
      return(1 - (X[3] - X[1])^2 / (X[3] - X[2]) / (X[3] - X[4]))
    } else if(X[1] >= X[3])
    {
      return(1)
    }
  }

  k <- max(q1, a1, b1, c1)
  if(k==1) return(pTest(c(lq, la, lb, lc)))

  params <- matrix(nrow=k, ncol=4)
  tryCatch(
  {
    params[,1] <- lq
    params[,2] <- la
    params[,3] <- lb
    params[,4] <- lc
  }, error = function(X) {
    stop(paste(" -- Argument Lengths: length of q = ", q1,
                ", a = ", a1, ", b = ", b1, ", c = ", c1, " -- ", X, sep=""))
  })

  return(apply(params, 1, pTest))
}

qltriangle <- function(p, a=1, b=100, c=10^((log10(a)+log10(b))/2), logbase=10)
{
  p1 <- length(p)
  a1 <- length(a)
  b1 <- length(b)
  c1 <- length(c)

  stopifnot(length(logbase) == 1)
  if (logbase==10)
  {
    la <- log10(a)
    lb <- log10(b)
    lc <- log10(c)
  } else
  {
    la <- log(a)/log(logbase)
    lb <- log(b)/log(logbase)
    lc <- log(c)/log(logbase)
  }

  qTest <- function(X)
  {
    # X = c(p, a, b, c)
    if(any(is.na(X))){ # is.na is TRUE for NA, NaN, and FALSE
      if(any(is.nan(X))) return(NaN) # to conform to qunif
      else return(NA) # to conform to qunif
    } else if(X[2] > X[4] | X[3] < X[4])
    {
      warning("values required to be  a <= c <= b (at least one strict inequality)")
      return(NaN) # to conform to behavior of qunif
    } else if(X[1] < 0 | X[1] > 1)
    {
      warning("at least one p is outside [0,1]")
      return(NaN) # to conform to behavior of qunif
    } else if(any(is.infinite(X)))
    {
      return(NaN)
    } else if((X[2] != X[4] &&
               (X[2] + sqrt(X[1]*(X[3]-X[2])*(X[4]-X[2]))) <= X[4]) |
              (X[2] == X[4] &&
               (X[2] + sqrt(X[1]*(X[3]-X[2])*(X[4]-X[2]))) < X[4]))
    {
      return(logbase^(X[2] + sqrt(X[1]*(X[3]-X[2])*(X[4]-X[2]))))
    } else if((X[2] != X[4] &&
               (X[3] - sqrt((1-X[1])*(X[3]-X[2])*(X[3]-X[4]))) > X[4]) |
              (X[2] == X[4] &&
               (X[3] - sqrt((1-X[1])*(X[3]-X[2])*(X[3]-X[4]))) >= X[4]))
    {
      return(logbase^(X[3] - sqrt((1-X[1])*(X[3]-X[2])*(X[3]-X[4]))))
    } else stop("Unexpected Result")
  }

  k <- max(p1, a1, b1, c1)
  if(k==1) return(qTest(c(p,la,lb,lc)))

  params <- matrix(nrow=k, ncol=4)
  tryCatch(
  {
    params[,1] <- p
    params[,2] <- la
    params[,3] <- lb
    params[,4] <- lc
  }, error = function(X) {
    stop(paste(" -- Argument Lengths: length of p = ", p1,
                ", a = ", a1, ", b = ", b1, ", c = ", c1, " -- ", X, sep=""))
  })

  return(apply(params, 1, qTest))
}

