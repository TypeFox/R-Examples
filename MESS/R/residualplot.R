residualplot.default <- function(x, y, candy=TRUE, bandwidth = 0.3, xlab="Fitted values", ylab="Std.res.", col.sd="blue", col.alpha=0.3,...) { 

  if (candy) 
    plot(x, y, xlab = xlab, ylab = ylab, 
         pch = 1 + 15 * (abs(y) > 1.96), ...)
  else plot(x, y, xlab = xlab, ylab=ylab, ...)
  if (candy) {

    # Set the colors
    if (col.alpha == 0) 
      col.trans <- col.sd
    else col.trans <- sapply(col.sd, FUN = function(x) do.call(rgb, 
                                     as.list(c(col2rgb(x)/255, col.alpha))))
    uniqx <- sort(unique(x))
    if (length(uniqx) > 3) {
      lines(smooth.spline(x, y, df = 3), lty = 2, lwd = 2, 
            col = "black")
    }
    window <- bandwidth * (max(x) - min(x))/2
    vary <- length(uniqx)
    for (i in 1:length(uniqx)) {
      vary[i] <- 1.96 * sd(y[abs(x - uniqx[i]) <= window])
    }
    vary[is.na(vary)] <- 0
    polygon(c(uniqx, rev(uniqx)), c(vary, -(rev(vary))), 
            col = col.trans, border = NA)
  }
  return(invisible(NULL))
}


residualplot.lm <- function(x, y, candy=TRUE, bandwidth = 0.3, xlab="Fitted values", ylab="Stud.res.", col.sd="blue", col.alpha=0.3,...) {
  y <- rstudent(x)
  x <- predict(x)
  residualplot(x, y, candy, bandwidth, xlab, ylab, col.sd, col.alpha, ...)
}



residualplot <- function(x, y, candy=TRUE, bandwidth = 0.3, 
	                 xlab="Fitted values", ylab="Std.res.", 
                         col.sd="blue", col.alpha=0.3,...) {
  UseMethod("residualplot")
}


