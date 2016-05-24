pselect <- function(n, p, min.diff=NULL, min.resp=NULL) {
  ntrt <- length(p)
# check that there are at least 2 treatments and p is a probability vector
  if (ntrt <= 1) stop("there should be at least 2 treatments")
  if (min(p) < 0 | max(p) > 1) stop("p should be a vector of probabilities i.e. between 0 and 1")
# unequal sample size is allowed for comparing two treatments only
  nlen <- length(n)
  if (nlen > 2) stop("sample size should be a vector of size 1 or 2")
  if (nlen==2 & n[1] == n[2]) {
    n <- n[1]
    nlen <- 1
  }
  out <- list()
  if (nlen == 1) {
    if (missing(min.resp)) min.resp <- 0
    if (missing(min.diff)) min.diff <- 1
    if (min.diff <= 0) stop("min.diff should be positive")
    if (min.diff < 1) min.diff <- ceiling(n*min.diff)
    if (min.diff != round(min.diff)) stop("if min.diff > 1 it should be a positive integer")
    psel0 <- prod(pbinom(min.resp-1, n, p))
    psel <- rep(0, ntrt)
    for(i in min.resp:n) {
      i0 <- max(i-min.diff, min.resp-1)
      if (i0 >= 0) {
        pb <- pbinom(i0, n, p)
        psel <- psel + dbinom(i,n,p)*prod(pb)/pb
      }
    }
    if (min.resp > 0)  out$prob.none.worthy <- psel0
  }
  if (nlen ==2) {
    if (ntrt > 2) stop("n and p are not of equal length")
    if (missing(min.diff)) stop ("min.diff should be specified as a rate for unequal sample size")
    if (min.diff <= 0 | min.diff >= 1) stop ("min.diff should be in (0,1) for unequal sample size")
    if (missing(min.resp)) {
      min.resp <- rep(0, 2)
    }
    if (length(min.resp) != 2) stop("min.resp should be the same length as n (one for each treatment)")
    if (any(min.resp != round(min.resp)) | min(min.resp) < 0 ) stop("min.resp should be a non-negative integer")
    n1 <- n[1]
    n2 <- n[2]
    p1 <- p[1]
    p2 <- p[2]
    psel0 <- 0
    psel <- rep(0, 2)
    for (i in 0:n1) {
      for (j in 0:n2) {
        pij <- prod(dbinom(c(i,j), n, p))
        if (i < min.resp[1] & j < min.resp[2]) {
          psel0 <- psel0 + pij
        } else {
          if ((i/n1 - j/n2 >= min.diff) | (j < min.resp[2])) psel[1] <- psel[1] + pij
          if ((j/n2 - i/n1 >= min.diff) | (i < min.resp[1])) psel[2] <- psel[2] + pij
        }
      }
    }
    if (max(min.resp) > 0)  out$prob.none.worthy <- psel0
  }
  out$prob.inconclusive <- 1-sum(psel)-psel0
  out$prob.selection <- cbind(n, p, psel)
  colnames(out$prob.selection) <- c("sample.size", "resp.rate", "prob.selection")
  rownames(out$prob.selection) <- paste("trt", 1:ntrt)
  out
}
