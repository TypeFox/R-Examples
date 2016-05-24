# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


freqresponse = function (b, a, fs = 0, add = FALSE, show = TRUE, steps = 1000,...){
  if (fs > 0) stepsize = (2/fs)
  w = seq (0, pi, length.out = steps)
  j = sqrt (as.complex(-1))
  numerator = 0
  denomenator = 0
  
  for (i in 1:length(b)) numerator = numerator + b[i] * (cos((i-1) * w) - sin((i-1) * w) * j)
  if (length(a) > 1)
  for (i in 2:length(a)) denomenator = denomenator - a[i] * (cos((i-1) * w) - sin((i-1) * w) * j)

  denomenator = 1 - denomenator
  
  response = log(abs (numerator / (denomenator)), 10) * 20
  response = response - max (response)
  
  if (add == FALSE & show == TRUE){
    if (fs == 0) plot (w / (pi*2), response, xlab = '(Frequency / Sampling Freq.)', ylab = 'Power (dB)',
                       type = 'l', lwd = 2, xaxs = 'i', ...)
    if (fs > 0) plot (w*((fs/2) / pi), response, xlab = 'Hz', ylab = 'Power (dB)',
                      type = 'l', lwd = 2, xaxs = 'i', ...)
  }
  if (add == TRUE & show == TRUE){
    if (fs == 0) lines (w / (pi*2), response, lwd = 2, ...)
    if (fs > 0) lines (w*((fs/2) / pi), response, lwd = 2, ...)
  }  
  if (fs == 0) out = data.frame (frequency = w / (pi*2), response = response)
  if (fs > 0) out = data.frame (frequency = w*((fs/2) / pi), response = response)
  invisible (out)
}
