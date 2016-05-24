`plot.biso` <- function(x, xnew, cl=TRUE, add=FALSE, color="red", ...){
  o <- order(x$x)
  z <- x$x[o]
  y <- x$y[o]
  m <- x$m
  u <- as.matrix(x$postdraws[,1:m])
  u0 <- x$postdraws$u0
  if (!add) plot(z, y, type="n", ...)
  if (missing(xnew)) xnew <- z
  wnew <- get.W((xnew - x$xmin)/x$xrng, m)
  if (cl){
    lcl <- sapply(1:nrow(wnew), function(i) quantile(u0 + u%*%wnew[i,], probs=0.025))
    ucl <- sapply(1:nrow(wnew), function(i) quantile(u0 + u%*%wnew[i,], probs=0.975))
    polygon( c(xnew,rev(xnew)), c(lcl,rev(ucl)), col=gray(0.75), border=NA)
    }
  points(z, y)
  yhat <- sapply(1:nrow(wnew), function(i) median(u0 + u%*%wnew[i,]))
  lines(xnew, yhat, col=color, ...)
  }

