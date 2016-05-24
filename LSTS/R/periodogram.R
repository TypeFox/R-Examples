periodogram = function(y, plot = TRUE, include.taper = FALSE){
  series = y - mean(y, na.rm = TRUE)
  N = sum(is.na(series))
  series[is.na(series)]=0
  n = length(series)
  if(include.taper == TRUE){
    a = 0.5*(1-cos(2*pi*(0:(n-1))/n))
    series = series*a
  }
  aux = Mod(fft(series))^2
  m = n/2
  periodogram = (aux[2:(m+1)])/(2*pi*(n-N))
  if(include.taper == TRUE){
    periodogram = (aux[2:(m+1)])/(3*pi*(n-N)/4)	
  }
  lambda = (2*pi * (1:m))/n
  if(plot == TRUE){
    plot(periodogram ~ lambda, bty = "n", las = 1, xlab = expression("Frequency"), ylab = expression("Periodogram"), xaxt = "n", type = "l")
    axis(1, at = seq(0,pi,pi/4), labels = expression(0, pi/4, pi/2, 3*pi/4, pi))
  }
  return(list(periodogram = periodogram, lambda = lambda))
}
