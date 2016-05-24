`plot.cardiPeakwindow` <-
function (x, y, add=FALSE, ...){
  peaks <- x 
  x <- peaks$data$x
  y <- peaks$data$y
  smd         <- peaks$smd.indices
  thecol      <- rep(3, length(peaks$peakid))
  thecol[smd] <- 2
  if (add) {
    lines(x[smd], y[smd], col = "red", lwd = 2)
    text(x[smd], y[smd], peaks$peakid[smd], col = "red", pos = 3)
  } else {
    plot(x, y, type = "l", ...)
    lines(x[smd], y[smd], col = "red", lwd = 2)
    text(x, y, peaks$peakid, col = thecol, pos = 3)
  }
}

