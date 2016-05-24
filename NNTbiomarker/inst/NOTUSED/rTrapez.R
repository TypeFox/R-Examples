rTrapez = function(n, slope) {
  #must be between -2 and 2.  heights run from 1-slope/2 to 1+slope/2
  #   F = (1-slope/2)*x + slope*x^2/2
  #   Finv = -(1-slope/2)+sqrt((1-slope/2)^2 - 4*slope/2*(-F))/(2*slope/2)
  if(is.infinite(1/slope)) slope = 1e-7
  return((-(1-slope/2)+sqrt((1-slope/2)^2 + 2*slope*runif(n)))/slope)
}
