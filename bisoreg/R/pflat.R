`pflat` <- function(obj,...) UseMethod("pflat")

`pflat.default` <- function(obj,...) stop("No default method for pflat.  Tough break.")

`pflat.biso` <- function(obj,...){
  M <- obj$m
  draws <- obj$postdraws
  mean(rowSums(draws[,1:M])==0)
  }

