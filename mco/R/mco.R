##
## mco.R - General MCO utilities
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.uni-dortmund.de>
##

## Extract pareto set:
paretoSet <- function(x, ...) 
  UseMethod("paretoSet", x)

## Extract pareto front:
paretoFront <- function(x, ...)
  UseMethod("paretoFront", x)

## Short helper functions to return matricies as is...
paretoSet.matrix <- function(x, ...)   { return(x) }
paretoFront.matrix <- function(x, ...) { return(x) }

paretoFilter <- function(x, ...)
  UseMethod("paretoFilter",x)

paretoFilter.matrix <- function(x, ...) {
  d <- ncol(x)
  n <- nrow(x)
  is.optimal <- rep(TRUE, n)
  for(i in 1:(n-1)) {
    for (j in i:n) {
      if (i != j && (is.optimal[i] || is.optimal[j])) {
        xi <- x[i,]
        xj <- x[j,]
        if (all(xi <= xj) && any(xi < xj)) { ## i dominates j
          is.optimal[j] <- FALSE
        } else if (all(xj <= xi) && any(xj < xi)) { ## j dominates i
          is.optimal[i] <- FALSE
        }
      }
    }
  }
  ## Always return a matrix!
  return(x[is.optimal,,drop=FALSE])
}
