# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


windowfunc = function (npoints, type = 'hann', parameter = -1){
  if (length (npoints) > 1) npoints = length (npoints) 

  if (!is.numeric (npoints)) stop ('Invalid number of points specified.')
  npoints = round (npoints)
  N = npoints
  n = 0:(N-1)
  output = NULL
  
  if (type == 'rectangular') output = rep (1, N)
  if (type == 'blackman') output = (7938/18608) - (9240/18608)*cos((2*pi*n)/(N-1)) + (1430/18608)*cos((4*pi*n)/(N-1))
  if (type == 'hann' | type == 'hanning') output = 0.5 * (1 - cos ((2*pi*n)/(N-1)))
  if (type == 'hamming') output = 0.54 - 0.46 * cos ((2*pi*n)/(N-1))
  if (type == 'cosine' | type == 'sine') output = sin ((n*pi) / (N-1))
  if (type == 'bartlett') output = (2 / (N-1)) * ((N-1)/2 - abs(n - (N-1)/2))
  if (type == 'gaussian'){
    if (parameter == -1) parameter = 0.4
    output = exp (-.5 * ((n - (N-1)/2) / (parameter * (N-1)/ 2))^2)
  }
  if (type == 'kaiser'){
    if (parameter == -1) parameter = 2
    output = besselI (parameter*pi * sqrt (1 - (2*(n)/(N-1) -1)^2), 0) / besselI(parameter*pi, 0)
  }
  if (type =='') stop ('No window type provided.')
  if (is.null(output)) stop ('Invalid window type provided.')
  return (output)
}

