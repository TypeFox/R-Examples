##
## minimum_di.R - minimum desirability index
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.tu-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.tu-dortmund.de>
##

minimumDI <- function(f, ...)
  UseMethod("minimumDI", f)

## Vector input
minimumDI.numeric <- function(f, ...) 
  min(f)

## Matrix input
minimumDI.matrix <- function(f, margin=1, ...)
  apply(f, margin, min)

## Array input
minimumDI.array <- function(f, margin=1, ...)
  apply(f, margin, min)

  
minimumDI.desire.function <- function(f, ...) {
  ev <- function(x) {
    fn <- function(z)
      min(sapply(i, function(k) dfs[[k]](z[k])))

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
  class(ev) <- c("desire.index", "mimimumDI")
  return(ev)
}

minimumDI.composite.desire.function <- function(f, ...) {
  ev <- function(x) {
    dval <- sapply(dfs, function(f) f(x))
    
    if (is.matrix(dval)) {
      apply(dval, 1, min)
    } else {
      min(dval)
    }
  }

  dfs <- list(f, ...)
  if (!all(sapply(dfs, is.composite.desirability)))
    stop("Not all supplied arguments are composite desirability functions.")
  
  class(ev) <- "desire.index"
  return(ev)
}
