"qqexp" <-
function(y, line=FALSE, ...)
{ y <- y[!is.na(y)]
  n <- length(y)
  x <- qexp(c(1:n)/(n+1))
  m <- mean(y)
  if (any(range(y)<0)) stop("Data contain negative values")
  ylim <- c(0,max(y))
  qqplot(x, y, xlab="Exponential plotting position",ylim=ylim,ylab="Ordered sample", ...)
  if (line) abline(0,m,lty=2)
  invisible()
}

