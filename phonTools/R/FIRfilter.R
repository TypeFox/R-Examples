# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


FIRfilter = function (sound, from = 0, to = fs/2, fs = 22050, order = 200, 
                      verify = FALSE, impulse = NULL, pad = TRUE){
  soundout = 0; tsout = 0;	
  if (order%%2) order = order +1
  if (class(sound) == "ts"){
    fs = frequency(sound)
    tsout = 1
  }
  if (class(sound) == "sound") {
    soundout = 1
    tmp = sound
    fs = sound$fs
    sound = sound$sound
  }
  if (order > length (sound)) order = length(sound)
  if (order < 4) order = 4
  if (from < 0) stop ('Low cutoff must be above 0 Hz.')
  if (from > fs/2) stop ('High cutoff must be below the Nyquist frequency.')
  sdin = sd (sound)
  
  maxamp = max(sound)
  M = order
  n = seq(0,M-1,1)
  Mi = (M-1) / 2  
  fromrads = (((fs/2)-from) / fs) * pi
  torads = (to / fs) * pi
  fromrads2 = (from / fs) * pi
  
  fromh = (-1)^(n)*2*fromrads*sinc((2*fromrads)*(n-Mi))  ##min freq passed
  toh = 2*torads*sinc(2*torads*(n-Mi))  ##max freq passed
  fromh2 = 2*fromrads2*sinc(2*fromrads2*(n-Mi))  
  
  fromh = fromh * windowfunc(length(fromh), type ='blackman')
  toh = toh * windowfunc(length(toh), type ='blackman')
  fromh2 = fromh2 * windowfunc(length(fromh2), type ='blackman')
  
  if (pad) sound = c(rep(sound[1], order/2), sound, rep(tail(sound,1),order/2))
  if (from!=0 & to==fs/2) output = filter (sound, fromh, method = 'convolution')
  if (from==0 & to!=fs/2) output = filter (sound, toh, method = 'convolution')
  if (from!=0 & to!=fs/2) output = filter (sound, fromh2-toh, method = 'convolution')
  if (!is.null(impulse)) output = filter (sound, impulse, method = 'convolution')
  output = output[!is.na(output)]
  output = output/sd(output) * sdin
  if (pad) output = output[-1]
  
  if (verify == TRUE){
    spectralslice (sound, fs = fs, col = 3, lty = 'dashed', ylim = c(-110,0), padding = 0, window = 'kaiser')  
    spectralslice (output, fs = fs, add = TRUE, padding = 0, window = 'kaiser') 
    abline (v = c(from,to), lwd = 2, col = 2)
  }  
  if (soundout == 1){
    tmp$sound = output
    output = tmp
  }
  else if (tsout == 1){
    output = ts (output, frequency = fs, start = 0)
  }
  invisible (output)
}
