

c.rv <- function(..., recursive=FALSE)
{
  ## a kludge to disable dispatching rv
  x <- c(list(NA), ..., recursive=recursive)[-1]
  class(x) <- class(rv())
  return(x)
}

c.rvsummary <- function(..., recursive=FALSE)
{
  ## a kludge to disable dispatching rv
  x <- c(list(NA), ..., recursive=recursive)[-1]
  class(x) <- "rvsummary"
  return(x)
}

cc <- function(..., recursive=FALSE) {
  args <- lapply(as.list(match.call())[-1], eval, envir=parent.frame())
  cls <- sapply(args, class)
  if ("rvsummary" %in% cls) {
    x <- c.rvsummary(..., recursive=recursive)
  } else if ("rv" %in% cls) {
    x <- c.rv(...)
  } else {
    x <- c(..., recursive=recursive)
  }
  return(x)
}



