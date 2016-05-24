
# FUNCTION:  KERNEL EPANECHNIKOV
Epanechnikov <- function(u)
{ 
  1 - (u)^2
}


# FUNCTION:  KERNEL EPANECHNIKOV
quadratic <- function(u)
{ 
  1 - (u)^2
}


# FUNCTION:  GAUSSIAN KERNEL
gaussian <- function(u)
{ 
  (1/sqrt(2*pi))*exp((-1/2)*u^2)
}


# FUNCTION:  UNIFORM KERNEL
uniform <- function(u)
{ 
  num_row <- nrow(u)
  unif <- function(a) {
    if ((a < 1) & (a > -1)) 1
    else 0
  }
  matrix(sapply(u, unif), nrow=num_row)
}


# FUNCTION:  TRIWEIGHT KERNEL
triweight <- function(u)
{ 
  (1 - (u)^2)^3
}

