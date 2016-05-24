
rvmean <- function (x) {
  UseMethod("rvmean")
}

rvmean.rv <- function (x) {
  m <- colMeans(sims(x), na.rm=TRUE)
  names(m) <- names(x)
  dim(m) <- dim(x)
  dimnames(m) <- dimnames(x)
  return(m)
}

rvmean.rvsummary <- function (x) {
  unlist(rvattr(x, "mean"), use.names=TRUE)
}

rvmean.default <- function (x) {
  if (!is.numeric(x)) {
    x[] <- as.numeric(x)
  }
  return(x)
}

E <- function (x) {
  rvmean(x)
}

Pr <- function (X) {
  if (! is.logical.rv(X)) {
      stop("Argument for Pr must be a logical statement such as 'x>0'")
  }
  rvmean(X)
}
