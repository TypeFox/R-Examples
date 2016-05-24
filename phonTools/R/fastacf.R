# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


fastacf = function (signal, lag.max = length(signal), window = 'hann',
                    show = TRUE, correct = FALSE){
  if (class(signal)=='sound') signal = signal$sound
 
  n = length (signal)
  if (lag.max > n) lag.max = n 
  
  signal = signal - mean (signal)
  signal = c(signal*windowfunc(signal, type = window), zeros (signal))
  
  pad=0
  #if (pad){  
  #  pad = 2^ceiling (log (length (signal),2)) - length (signal)
  #  signal = c(signal, zeros (pad))
  #}
  s = fft (signal)
  s = Re (fft (s * Conj(s), inverse = TRUE))
  s = s[1:(lag.max)] / s[1]
  
  if (correct){
    w = fft (c(windowfunc (n, type = window), zeros(n), zeros(pad)))
    w = Re (fft (w * Conj(w), inverse = TRUE))
    w = w[1:(lag.max)] / w[1]
    s = s/w
  }
  lags = (0:(lag.max-1))
  
  if (show){
    plot (lags, s, ylim = c(-1,1), xlab = 'Lag', ylab = 'ACF', type = 'l',xaxs = 'i')
    if (correct) abline (v = n - 50, col = 2)
    abline (h = 0, lty = 'dotted')
  }  
  invisible (data.frame (lag = lags, acf = s))
}

