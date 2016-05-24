# Copied directly from seewave 1.7.3
# Modified: 2014 MAR 20

hanning.w <- function (n)
{
  if(n <= 0) stop("'n' must be a positive integer")
  n <- n-1
  w <- 0.5-0.5*cos(2*pi*(0:n)/n)
  return(w)
}


hamming.w <- function (n)
{
  if(n <= 0) stop("'n' must be a positive integer")
  n <- n-1
  w <- 0.54-0.46*cos(2*pi*(0:n)/n)
  return(w)
}

blackman.w <- function (n)
{
  if(n <= 0) stop("'n' must be a positive integer")
  n <- n-1
  w <- 0.42-0.5*cos(2*pi*(0:n)/n)+0.08*cos(4*pi*(0:n)/n)
  return(w)
}

flattop.w <- function (n)
{
  if(n <= 0) stop("'n' must be a positive integer")
  n <- n-1
  w <- 0.2156-0.4160*cos(2*pi*(0:n)/n)+0.2781*cos(4*pi*(0:n)/n)
  -0.0836*cos(6*pi*(0:n)/n)+0.0069*cos(8*pi*(0:n)/n)   
  return(w)
}

rectangle.w <- function (n)
{
  if(n <= 0) stop("'n' must be a positive integer")
  w <- rep(1, n)
  return(w)
}

bartlett.w <- function (n)
{
  if(n <= 0) stop("'n' must be a positive integer")

  n <- n-1
  m <- n%/%2
  w <- c((2*(0:(m-1)))/n, 2-((2*(m:n))/n))
  return(w)
}







