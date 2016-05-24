# Test functions ----------------------------------------------------------

# method for plotting (2D and 3D)
# method for summary and print
# f() print the minimum

summary.calibrar.function = function(object, ...) {
  return(invisible())
}

#' Sphere function with random noise 
#' 
#' This function calculates the Euclidian distance from a point to the origin
#' after a random displacement of it position. 
#' 
#' @param x The coordinates of the point
#' @param sd The standard deviation of the noise
#' to be added to the position of \code{x}, a normal distribution with mean
#' zero is used.
#' @param aggregate If \code{aggregate} is \code{TRUE} the distance is returned, 
#' otherwise the size of the projection of the distance among each axis.
#' @return The distance from the point \code{x} to the origin after a random
#' displacement.
#' @author Ricardo Oliveros--Ramos
#' @keywords stochastic random
#' @examples
#' SphereN(rep(0, 10))
#' @export
SphereN = function(x, sd=0.1, aggregate=TRUE) {
  # f(0,...,0) = 0
  # x_i \in ]-Inf, Inf[
  x = x + rnorm(length(x), sd=sd)
  out = x^2
  if(isTRUE(aggregate)) return(sum(out)) else return(out) 
}
