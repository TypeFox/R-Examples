w.median <- function(x,w) {
  if (missing(w)) w <- rep(1,length(x))
  ok <- complete.cases(x,w)
  x <- x[ok]
  w <- w[ok]
  ind <- sort.list(x)
  x <- x[ind]
  w <- w[ind]
  ind1 <- min(which(cumsum(w)/sum(w) >=0.5))
  ind2 <- if( (w[1] / sum(w)) > 0.5 ) { 1 } else {max(which(cumsum(w)/sum(w) <=0.5))}
  max(x[ind1],x[ind2])
}
