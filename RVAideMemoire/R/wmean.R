wmean <- function(x,w=rep(1,length(x)),na.rm=TRUE) {
  if (na.rm) {
    comp <- complete.cases(x,w)
    x <- x[comp]
    w <- w[comp]
  }
  sum(x*w)/sum(w)
}
