#' Bernstein polynomials
#' 
#' Bernstein polynomials
#' 
#' @param k Bernstein polynomial k
#' @param n Bernstein polynomial degree
#' @param indeterminate indeterminate
#' @return a mpoly object
#' @author David Kahle
#' @export
#' @examples
#' 
#' bernstein(0, 0)
#' 
#' bernstein(0, 1)
#' bernstein(1, 1)
#' 
#' bernstein(0, 1, "t")
#' 
#' bernstein(0:2, 2)
#' bernstein(0:3, 3)
#' bernstein(0:3, 3, "t")
#' 
#' 
#' bernstein(0:4, 4)
#' bernstein(0:10, 10)
#' bernstein(0:10, 10, "t")
#' bernstein(0:20, 20, "t")
#' 
#' \dontrun{  # visualize the bernstein polynomials
#' 
#' library(ggplot2); theme_set(theme_classic())
#' library(tidyr)
#' 
#' s <- seq(0, 1, length.out = 101)
#' N <- 10 # number of bernstein polynomials to plot
#' (bernPolys <- bernstein(0:N, N))
#' 
#' df <- data.frame(s, as.function(bernPolys)(s))
#' names(df) <- c("x", paste0("B_", 0:N))
#' head(df)
#' 
#' mdf <- gather(df, degree, value, -x)
#' head(mdf)
#' 
#' qplot(x, value, data = mdf, geom = "line", color = degree)
#' 
#' }
#' 
bernstein <- function(k, n, indeterminate = "x"){  
  
  ## make it possible for vector k args
  if(length(k) > 1){
    listOPolys <- lapply(k, function(.) bernstein(., n, indeterminate))
    class(listOPolys) <- "mpolyList"
    return(listOPolys)
  }
  
  ## construct coefficients and degrees of terms
  m <- n - k  
  coefs <- choose(n, k) * (-1)^(0:m) * choose(m, 0:m)
  degs  <- k:n
  
  ## construct polynomial as list
  p <- Map(function(deg, coef) c(x = deg, coef = coef), degs, coefs)
  
  ## wipe out zeros
  p <- lapply(p, function(v) v[v != 0])
  
  ## class list
  class(p) <- c("bernstein", "mpoly")
  attr(p, "bernstein") <- list(k = k, n = n, indeterminate = indeterminate)

  ## swap and return
  swap(p, "x", indeterminate)
}





















#' Bernstein polynomial approximation
#' 
#' Bernstein polynomial approximation
#' 
#' @param f the function to approximate
#' @param n Bernstein polynomial degree
#' @param lower lower bound for approximation
#' @param upper upper bound for approximation
#' @param indeterminate indeterminate
#' @return a mpoly object
#' @author David Kahle 
#' @export
#' @examples
#' 
#' 
#' 
#' 
#' 
#' \dontrun{  # visualize the bernstein polynomials
#' 
#' library(ggplot2); theme_set(theme_bw())
#' library(reshape2)
#' 
#' 
#' 
#' 
#' f <- function(x) sin(2*pi*x)
#' p <- bernsteinApprox(f, 20) 
#' round(p, 3)
#' 
#' x <- seq(0, 1, length.out = 101)
#' df <- data.frame(
#'   x = rep(x, 2), 
#'   y = c(f(x), as.function(p)(x)), 
#'   which = rep(c("actual", "approx"), each = 101)
#' )
#' qplot(x, y, data = df, geom = "line", color = which)
#' 
#' 
#' 
#' 
#' 
#' 
#' p <- bernsteinApprox(sin, 20, pi/2, 1.5*pi) 
#' round(p, 4)
#' 
#' x <- seq(0, 2*pi, length.out = 101)
#' df <- data.frame(
#'   x = rep(x, 2), 
#'   y = c(sin(x), as.function(p)(x)), 
#'   which = rep(c("actual", "approx"), each = 101)
#' )
#' qplot(x, y, data = df, geom = "line", color = which)
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' p <- bernsteinApprox(dnorm, 15, -1.25, 1.25) 
#' round(p, 4)
#' 
#' x <- seq(-3, 3, length.out = 101)
#' df <- data.frame(
#'   x = rep(x, 2), 
#'   y = c(dnorm(x), as.function(p)(x)), 
#'   which = rep(c("actual", "approx"), each = 101)
#' )
#' qplot(x, y, data = df, geom = "line", color = which)
#' 
#' 
#' 
#' 
#' 
#' 
#' }
#'
bernsteinApprox <- function(f, n, lower = 0, upper = 1, indeterminate = "x"){  

  ## compute support and determine weights
  s <- (0:n)/n
  fscaled <- function(.) f( (upper-lower)*. + lower )
  weights <- as.list(fscaled(s))
  
  ## convert weights to mpolyList
  weights <- lapply(weights, function(x) mpoly(list(c(coef = x))))
  class(weights) <- "mpolyList"
  
  ## multiply weights by basis
  approxPoly <- Reduce(`+`, weights * bernstein(0:n, n, "temp"))  
  
  ## compute plugin and plug in
  pluginPoly <- (upper-lower)^-1 * (mp(indeterminate) + -1*lower)
  plug(approxPoly, "temp", pluginPoly)
  
}






