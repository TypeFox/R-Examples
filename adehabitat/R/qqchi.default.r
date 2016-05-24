"qqchi.default" <-
function (y, df = 2, ylim, main = "Chi Q-Q Plot",
                           xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
                           plot.it = TRUE, datax = FALSE,...) 
{
  if (has.na <- any(ina <- is.na(y))) {
    yN <- y
    y <- y[!ina]
  }
  if (0 == (n <- length(y))) 
    stop("y is empty or has only NAs")
  if (plot.it && missing(ylim)) 
    ylim <- range(y)
  x <- qchi(ppoints(n),df)[order(order(y))]
  if (has.na) {
    y <- x
    x <- yN
    x[!ina] <- y
    y <- yN
  }
  if (plot.it) 
    if (datax) 
      plot(y, x, main = main, xlab = ylab, ylab = xlab, 
           xlim = ylim, ...)
    else plot(x, y, main = main, xlab = xlab, ylab = ylab, 
              ylim = ylim, ...)
  invisible(if (datax) list(x = y, y = x) else list(x = x, 
                                     y = y))
}

