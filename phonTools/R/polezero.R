# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


polezero = function (b, a, ...){
  if (!is.numeric(a)) stop ('Inappropriate feedback coefficients.')
  if (!is.numeric(a)) stop ('Inappropriate feedforward coefficients.')

  aroots = polyroot (rev(a))
  broots = polyroot (rev(b))
  if (length(a) == 1) aroots = complex (real = 0, imaginary = 0)
  if (length(b) == 1) broots = complex (real = 0, imaginary = 0)

  xrange = yrange = range (-1.1,1.1)

  realmax = max(abs (Re(aroots)), abs (Re(broots)))
  imagmax = max(abs (Im(aroots)), abs (Im(broots)))
 
  if (realmax > 1) xrange = c(-realmax, realmax)
  if (imagmax > 1) yrange = c(-imagmax, imagmax)
 
  plot (aroots, xlim = xrange, ylim = yrange, pch = 4, lwd = 2, 
        xlab = 'Real', ylab = 'Imaginary', cex = 1.5)
  points (broots, lwd = 2, cex = 1.75)

  sdellipse (means = c(0,0), points = matrix (c(1,0,0,1),2,2), stdev = 1, density = .01)
  abline (h = 0, v = 0, lty = 'dotted')
}

