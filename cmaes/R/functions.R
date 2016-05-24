##
## functions.R - Test functions
##
## Author:
##  Olaf Mersmann (OME) <olafm@@statistik.tu-dortmund.de>
##

##'
##' Returns a new function
##' \deqn{g(x) = f(x - offset).}
##'
##' @param f test function
##' @param offset offset.
##' @return The shifted test function.
##' @export
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
shift_function <- function(f, offset) {
  force(offset);
  newf <- function(x)
    f(x - offset)
  return(newf)
}

##' Create a rotated test function
##'
##' Returns a new rotated test function defined as
##' \deqn{g(x) = f(Mx).}
##'
##' @param f test function.
##' @param M orthogonal square matrix defining the rotation.
##' @return The rotated test function.
##' @export
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
rotate_function <- function(f, M) {
  M <- t(M)
  newf <- function(x)
    f(drop(M %*% x))
  return(newf)
}

##' Create a biased test function
##'
##' Returns a new biased test function defined as
##' \deqn{g(x) = f(x) + bias.}
##'
##' @param f test function
##' @param bias bias value.
##' @return The biased test function.
##' @export
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
bias_function <- function(f, bias) {
  force(bias);
  newf <- function(x)
    f(x) + bias
  return(newf)
}

################################################################################
##' Sphere function
##'
##' \deqn{f(x) = x'x}
##'
##' @param x parameter vector.
##' @export
f_sphere <- function(x)
  crossprod(x)

################################################################################
##' Random function
##'
##' \deqn{f(x) = runif(1)}
##'
##' @param x parameter vector.
##' @export
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
f_rand <- function(x)
  runif(1)
 
##' Rosenbrock function
##'
##' @param x parameter vector.
##' @export
##' @author David Arnu \email{david.arnu@@tu-dortmund.de}
f_rosenbrock <- function(x) {
  d <- length(x)
  z <- x + 1
  hz <- z[1:(d-1)]
  tz <- z[2:d]
  s <- sum(100 * (hz^2 - tz)^2 + (hz - 1)^2)
  return(s)
}

################################################################################
##' Rastrigin function
##'
##' @param x parameter vector.
##' @export
##' @author David Arnu \email{david.arnu@@tu-dortmund.de}
f_rastrigin <- function(x)
  sum(x*x - 10 * cos(2*pi*x) + 10)
