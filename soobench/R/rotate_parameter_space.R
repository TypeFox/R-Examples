##' Generate a random \code{d}-dimensional rotation matrix.
##'
##' The algorithm used to randomly create the rotation matrix is due
##' to R Salomon (see reference). No guarantee is given that the
##' generated rotation matricies are uniformly distributed in any
##' sense.
##' 
##' @param d Dimension of desired rotation matrix.
##' @return A random \eqn{d \times d} rotation matrix.
##' @references
##' Salomon R. Re-evaluating genetic algorithm performance under coordinate
##' rotation of benchmark functions. A survey of some theoretical and practical
##' aspects of genetic algorithms. Biosystems. 1996;39(3):263-78. 
random_rotation_matrix <- function(d) {
  simple_rotation_matrix <- function(d, i, j, alpha) {
    R <- diag(d)
    R[i, i] <- cos(alpha)
    R[i, j] <- sin(alpha)
    R[j, i] <- -sin(alpha)
    R[j, j] <- cos(alpha)
    R
  }

  R <- diag(d)
  for (i in 2:d)
    R <- R %*% simple_rotation_matrix(d, 1, i, runif(1, -pi/4, pi/4))

  if (d > 2) {
    for (i in 2:(d-1))
      R <- R %*% simple_rotation_matrix(d, i, d, runif(1, -pi/4, pi/4))
  }
  R    
}

##' Rotate the parameter space of a SOO function.
##'
##' This function is a simple parameter space transformation. Given a
##' function \eqn{f(x)} it retuns a new function \eqn{f_r(x) = f(Rx)},
##' where \eqn{R} is a random rotation matrix.
##'
##' If you want repeatable results, make sure you explicitly set a
##' seed before calling \code{rotate_parameter_space}.
##' 
##' @param fn A \code{soo_function} object.
##' @return A new \code{soo_function} object where the parameter space
##'  has been randomly rotated.
##' @examples
##' f <- ackley_function(2)
##' f_r <- rotate_parameter_space(f)
##' par(mfrow=c(1, 2))
##' plot(f)
##' plot(f_r)
##' 
##' @export
rotate_parameter_space <- function(fn) {
  d <- number_of_parameters(fn)
  R <- random_rotation_matrix(d)
  opt <- global_minimum(fn)
  
  ## Handle multiple optima locations:
  opt$par <- if (is.list(opt$par)) {
    tmp <- Filter(function(x) is_in_bounds(fn, x),
                  lapply(opt$par, function(x) R %*% x))
    if (length(tmp) != length(opt$par)) {
      if (length(tmp) > 0)
        warning("Some global minima are not inside bounds after rotation.")
      else
        warning("No global minima are inside bounds after rotation.")
    }
    tmp
  } else {
    tmp <- R %*% opt$par
    if (!is_in_bounds(fn, tmp))
      warning("Global minimum is not inside bounds after rotation.")
    tmp
  }

  soo_function(name=sprintf("Rotated %s", attr(fn, "name")),
               id=sprintf("rotated-%s", function_id(fn)),
               fun=function(x, ...) fn(R %*% x),
               dimensions=d,
               lower_bounds=lower_bounds(fn),
               upper_bounds=upper_bounds(fn),
               best_par=opt$par,
               best_value=opt$value)
}

is_in_bounds <- function(fn, par) {
  lb <- lower_bounds(fn)
  ub <- upper_bounds(fn)
  all(lb <= par & par <= ub)
}
