plot.kerndwd = function(x, color=FALSE, ...) {
  alpha = x$alpha[-1, ]
  lambda = x$lambda
  index = log(lambda)
  iname = "log lambda"
  xlab = iname
  ylab = "coefficients"
  dotlist = list(...)
  type = dotlist$type
  if (is.null(type)) {
    if (color == FALSE) 
      matplot(index, t(alpha), lty=1, xlab=xlab, ylab=ylab, 
        type="l", pch=500, col=gray.colors(12, 
        start=0.05, end=0.7, gamma=2.2), ...) else 
      matplot(index, t(alpha), lty=1, xlab=xlab, ylab=ylab, 
        type="l", pch=500, ...)
  } else matplot(index, t(alpha), lty=1, xlab=xlab, ylab=ylab, ...)
} 
