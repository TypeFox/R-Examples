# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


interpolate = function (y, x = 1:length(y), steps = 20, increment = -1, show = FALSE, type = 'cubic', ...){
  if (length (y) < 3) stop ('Interpolation requires at least three points.')
  if (!is.numeric (y) | !is.numeric (x)) stop ('Non-numeric arguments provided')
  if (type == 'cubic'){
    n = length(y)
    h = x[2:n] - x[1:(n - 1)]
    A = matrix(0, n, n)
    A[1, ] = c(1, rep(0, n - 1))
    A[n, ] = c(rep(0, n - 1), 1)
    for (i in 2:(n - 1)) {
      A[i, i - 1] = h[i - 1]
      A[i, i] = 2 * (h[i - 1] + h[i])
      A[i, i + 1] = h[i]
    }
    r = matrix(0, n, 1)
    for (i in 2:(n - 1)) r[i] = 6 * (((y[i + 1] - y[i])/h[i]) - 
      ((y[i] - y[i - 1])/h[i - 1]))
    m = c(0, lm(r ~ A)$coefficients[3:n], 0)
    yinterp = NULL
    xinterp = NULL
	
    for (i in 1:(length(x) - 1)) {
      a = y[i]
      b = ((y[i + 1] - y[i])/h[i]) - ((h[i]/2) * m[i]) - (h[i]/6) * (m[i + 1] - m[i])
      c = m[i]/2
      d = (m[i + 1] - m[i])/(6 * h[i])
      if (increment <= 0) xx = seq(x[i], x[i + 1], length.out = steps)
      if (increment > 0) xx = seq(x[i], x[i + 1], by = increment)
      yy = a + b * (xx - x[i]) + c * ((xx - x[i])^2) + d * ((xx - x[i])^3)
      xinterp = c(xinterp, xx)
      yinterp = c(yinterp, yy)
    }
  }
  if (type == 'linear'){
    n = length (y) - 1
    xinterp = NULL
    yinterp = NULL
    for (i in 1:n){
      if (increment <= 0){
        xinterp = c(xinterp, seq (x[i], x[i+1], length.out = steps))
        yinterp = c(yinterp, seq (y[i], y[i+1], length.out = steps))
      }
      if (increment > 0){
        xinterp = c(xinterp, seq (x[i], x[i+1], by = increment))
        yinterp = c(yinterp, seq (y[i], y[i+1], length.out = length(seq (x[i], x[i+1], by = increment))))
        
      }
    }
  }
  if (show == TRUE) {
    plot(xinterp, yinterp, type = "l")
    points(x, y)
  }
  output = data.frame(x = xinterp, y = yinterp)
  output = output[output[,1] != c(output[-1,1], 'q'),]
  rownames(output) = 1:nrow(output)
  invisible (output)
}
