
obslength <- function (x) {
  UseMethod("obslength")
}

obslength.rv <- function (x) {
  s <- (! is.na(sims(x)))
  a <- apply(s, 1, which)
  sapply(a, function (x) {
    if (length(x) == 0) { 0 } else { max(x) }
  })
}

obslength.numeric <- function (x) {
  w <- which(! is.na(x))
  if (length(w) == 0) {
    return(0)
  } else {
    return(max(w))
  }
}
