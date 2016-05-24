# Some tests.  The functions are defined first, then executed on the
# last few lines.  The functions all have obvious names.

require(contfrac)

TOLERANCE <- 1e-10

sqrt.of.11 <- function(...){
  lhs <- sqrt(11)
  rhs <- CF(c(3,rep(c(3,6),30)))
  return(lhs - rhs)
}

"sqrt.of.71" <- function(...){
  lhs <- sqrt(71)
  rhs <- CF(c(8,rep(c(2,2,1,7,1,2,2,16),10)))
  return(lhs-rhs)
}

"exp.1" <- function(...){
  lhs <- exp(1)
  jj <- t(cbind(1,(1:30)*2,1))
  rhs <- CF(c(2,jj))
  return(rhs-lhs)
}

"eulergamma" <- function(...){
  lhs <- -psigamma(1)
  rhs <- CF(c(
              0, 1, 1, 2, 1, 2, 1, 4, 3, 13, 5, 1, 1, 8, 1, 2, 4, 1, 1, 40, 1, 11,
              3, 7, 1, 7, 1, 1, 5, 1, 49, 4, 1, 65, 1, 4, 7, 11, 1, 399, 2, 1, 3, 2,
              1, 2, 1, 5, 3, 2, 1, 10, 1, 1, 1, 1, 2, 1, 1, 3, 1, 4, 1, 1, 2, 5, 1,
              3, 6, 2, 1, 2, 1, 1, 1, 2, 1, 3, 16, 8, 1, 1, 2, 16
              ))  # CF expansion courtesy of Maple
  return(lhs-rhs)
}

"tan.1.1i" <- function(...){
  lhs <- tan(1 + 1i)
  z <- 1 + 1i
  rhs <- GCF(c(z,rep(-z^2,99)),seq(from=1,by=2,len=100))
  return(rhs-lhs)
}


stopifnot(all(as_cf(sqrt(2))[-1] == 2))

test_vector <- abs(c(sqrt.of.11() , sqrt.of.71() , exp.1() , eulergamma() , tan.1.1i() ))
stopifnot(all(test_vector < TOLERANCE))
