# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

makesound = function (sound, filename, fs = 22050){
  if (missing(filename))filename = paste (deparse(substitute(sound)), '.wav', sep='')
  if (!is.numeric(sound)) stop('The sound must be a numeric vector.')

  numSamples = length(sound)
  output = list(filename = filename, fs = fs, numSamples = numSamples, 
  duration = numSamples/fs * 1000, sound = ts(sound, frequency = fs, start=0))
  class(output) = "sound"
  output
}
