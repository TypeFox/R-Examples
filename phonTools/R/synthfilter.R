# Copyright (c) 2014 Santiago Barreda
# All rights reserved.


synthfilter = function (sound, band = c(0,fs/4), fs = 1, verify = FALSE, attenuation = 0){
  soundout = 0; tsout = 0;	
  if (class(sound) == "ts"){
    fs = frequency(sound)
    tsout = 1
  }
  if (class(sound) == "sound") {
    soundout = 1
    fs = sound$fs
    oldsound = sound
    sound = sound$sound
  } 
  sdin = sd (sound)
  n = length (sound)
  if (attenuation != 0) attenuation = 10^(-abs(attenuation)/10) * 2
  
  freqs = seq (attenuation, fs, length.out = n+1)[-(n+1)]
  pass = rep (attenuation, length(freqs))
  pass[freqs >= band[1] & freqs <= band[2]] = 1
  pass[freqs <= (fs-band[1]) & freqs >= (fs-band[2])] = 1
    
  soundspect = fft (sound)
  output = soundspect * pass                 
  output = Re(fft(output, inverse = TRUE))
  
  if (verify){
    par (mfrow = c(2,1), mar = c(4.5,4.5,3,1))
    spectralslice (sound, fs = fs, padding = 1000, main = 'Input Signal')
    spectralslice (output, fs = fs, padding = 1000, main = 'Output Signal')
  }
  output = output/sd(output) * sdin
  if (soundout == 1)  output = makesound (output, filename = oldsound$filename, fs = fs)
  if (tsout == 1) output = ts (output, frequency = fs, start = 0)
 
  invisible (output)
}

