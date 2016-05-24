plot.gpava <- function (x, main = "PAVA Plot", xlab = "Predictor", ylab = "Response", 
col = "lightblue",...)
{
#x ... object of class pava

  o <- order(x$z)
  xval <- x$z[o]
  yval <- x$x[o]
  xcum <-  c(xval[1] - mean(diff(xval)), xval) 
  jumps <- ((1:length(yval))[!duplicated(yval)]-1)[-1]   #jumps of fitted step function
  jumps <- c(1, jumps, length(xval))
  
  if (is.list(x$y)) {
    if (is.na(x$p)) {
      plot(xcum, c(NA, (mapply(x$solver, x$y, x$w)[o])), xlab = xlab, ylab = ylab, main = main)
    } else {
      plot(xcum, c(NA, (mapply(x$solver, x$y, x$w, x$p)[o])), xlab = xlab, ylab = ylab, main = main)
    }
  } else {
    plot(xcum, c(NA, x$y[o]), xlab = xlab, ylab = ylab, main = main)   
  }  
  lines(xval, yval, col = col, lwd = 1, type = "S")
  points(xval[jumps], yval[jumps], col = col, pch = 13)
  grid()
}
