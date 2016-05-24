residualplot <- function(object, bandwidth=.3, ...) {
  # Check that input is correct
  if (!inherits(object, "lm")) 
    stop("use only with \"lm\" objects")

  x <- predict(object)
  y <- rstandard(object)

  cc <- complete.cases(x, y)
  x <- x[cc]
  y <- y[cc]
  
  plot(x, y, 
       xlab="Fitted values", ylab="Standardized residuals", ...)
  # Make "outliers" black
  outliers <- (abs(y)>1.96)
  points(x[outliers], y[outliers], pch=16)

  # Add mean smoothing spline
  # If there is at least 4 unique x values
  uniqx <- sort(unique(x))
  if (length(uniqx)>3) {
    lines(smooth.spline(x, y, df=3), lty=2, lwd=2, col = "black")
  }

  # Slow approach here. Should be in c
  window <- bandwidth*(max(x)-min(x))/2
  vary <- length(uniqx)
  i <- 1
  for (j in uniqx) {
    vary[i] <- 1.96*sd(y[abs(x-j)<=window])
    i <- i +1
  }
  vary[is.na(vary)] <- 0
  color <- rgb(237, 149, 100, 50, maxColorValue=255)
  polygon(c(uniqx, rev(uniqx)), c(vary, -(rev(vary))),
          col=color, border=NA)
#  lines(uniqx, vary, col="red")

}
