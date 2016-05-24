smooth.periodogram = function(y, plot = TRUE, spar = 0){
  aux = periodogram(y, plot = FALSE)
  smooth.periodogram = smooth.spline(aux$periodogram, spar = spar)$y
  lambda = aux$lambda
  if(plot == TRUE){
    plot(smooth.periodogram ~ lambda, bty = "n", las = 1, xlab = expression("Frequency"), ylab = 	expression("Smooth Periodogram"), xaxt = "n", type = "l")
    axis(1, at = seq(0,pi,pi/4), labels = expression(0, pi/4, pi/2, 3*pi/4, pi))
  }
  return(list(smooth.periodogram = smooth.periodogram, lambda = lambda))
}
