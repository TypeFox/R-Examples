# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


preemphasis = function (input, cutoff = 50, fs = 22050, verify = FALSE, coeff = 0){
  soundout = 0
  tsout = 0
  if (class(input) == "ts"){
    fs = frequency(input)
    tsout = 1
  }
  if (class(input) == "sound") {
    soundout = 1
    oldsound = input
    fs = input$fs
    input = input$sound
  }
  if (coeff == 0) coeff = -exp (-2 * pi * cutoff / fs)
  out = as.numeric (filter (input, c(1, coeff), method = 'convolution', sides = 1))
  out[1] = input[1]
  if (verify == TRUE){
    spectralslice (out, fs = fs)
    spectralslice (input, fs = fs, add = TRUE, col = 3, lty = 'dotted')    
  }
  if (soundout == 1){
    oldsound$sound = out 
    invisible (oldsound)
  }
  else if (tsout == 1){
    out = ts (out, frequency = fs, start = 0)
    invisible (out)
  }
  else invisible (out)
}