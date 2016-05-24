# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


plot.rcr = function (x, y, arrangement = NULL, ...){
  coeffs = x$coefficients
  cols = ncol (coeffs) 
  if (is.null (arrangement)) arrangement = c(ceiling(cols/ 5), 5)
  
  oldpar = par()
  par (mfrow = arrangement, mar = c(4,4,2,1))
  for (i in 1:cols){
    mu = mean (coeffs[,i])
    stdev = sd (coeffs[,i])
    dense = density(coeffs[,i], kernel = 'gaussian')
    xrange = c(mu - 4*stdev, mu + 4*stdev)  
    xs = seq (xrange[1]-4*stdev, xrange[2]+4*stdev, abs(xrange[1]-xrange[2])/1000)
    ys = dnorm (xs, mu, stdev)
    yrange = c(0, max(ys, dense$y)*1.1)
    plot(dense$x, dense$y, main = colnames (coeffs)[i], xlim = xrange, ylim = yrange, type = 'l', 
    xlab = 'Coefficient', ylab = 'Density', cex.lab = 1.2, cex.axis = 1.2, lwd = 2, ...)
    abline (h = 0, v = 0)
    lines (xs, ys, col = 2, lty = 'dashed', lwd = 2)
  }  
  suppressWarnings (par (oldpar))
}
