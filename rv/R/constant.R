

# Value:  a logical vector (not rv), TRUE if a component is constant w.p. 1
#

is.constant <- function(x) {
  # Note: this corresponds to "constant with probability 1", while
  # rvvar(x)==0 would correspond to "constant almost surely"
  return(rvnsims(x)==1)
}

as.constant <- function(x)
{
  UseMethod('as.constant')
}

as.constant.rv <- function (x)
{
  z <- rvmean(x)
  z[! is.constant(x)] <- NA
  return(z)
}

as.constant.rvsummary <- function(x)
{
  for (i in which(is.constant(x))) {
    x[[i]] <- x[[i]][1]
  }
  for (i in which(!is.constant(x))) {
    x[[i]] <- NA
  }
  structure(unlist(x), names=names(x))  
}

as.constant.default <- function (x)
{
  return(x)
}

