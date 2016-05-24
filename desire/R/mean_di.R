##
## mean_di.R - weighted mean desirability index
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.tu-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.tu-dortmund.de>
##

meanDI <- function(f, ..., weights=1) 
  UseMethod("meanDI", f)

## Vector input
meanDI.numeric <- function(f, ..., weights) {
  weights <- weights/sum(weights)
  mean(f*weights)
}
  
## Matrix input
meanDI.matrix <- function(f, margin=1, ..., weights=1) {
  weights <- weights/sum(weights)
  apply(f, margin, function(x) mean(x*weights))
}

## Array input
meanDI.array <- function(f, margin=1, ..., weights=1)  {
  weights <- weights/sum(weights)
  apply(f, margin, function(x) mean(x * weights))
}

meanDI.desire.function <- function(f, ..., weights=1) {
  weights <- weights/sum(weights)
  ev <- function(x) {
    fn <- function(z)
      mean(sapply(i, function(k) dfs[[k]](z[k])) * weights)
    if (is.matrix(x)) {
      apply(x, 1, fn)
    } else {
      fn(x)
    }
  }
  
  dfs <- list(f, ...)
  if (!all(sapply(dfs, is.desirability)))
    stop("Not all supplied arguments are desirability functions.")
  
  i <- 1:length(dfs)
  class(ev) <- "desire.index"
  return(ev)
}

meanDI.composite.desire.function <- function(f, ..., weights=1) {  
  ev <- function(x) {
    dval <- sapply(dfs, function(f) f(x))

    if (is.matrix(dval)) {
      dval %*% weights
    } else {
      drop(crossprod(dval, weights))
    }
  }

  dfs <- list(f, ...)
  if (!all(sapply(dfs, is.composite.desirability)))
    stop("Not all supplied arguments are composite desirability functions.")

  l <- length(dfs)
  if (length(weights) == 1) 
    weights <- rep(weights, l)

  if (length(weights) != l)
    stop("Incorrect size for weights vector.")
  
  weights <- weights / sum(weights)  
  class(ev) <- "desire.index"
  return(ev)
}
