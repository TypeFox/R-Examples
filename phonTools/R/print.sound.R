# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

print.sound <-
function (x, ...){
  cat ("\n      Sound Object\n")
  cat ("\n   Read from file:        ", x$filename)
  cat ("\n   Sampling frequency:    ", x$fs, ' Hz')
  cat ("\n   Duration:              ", x$duration,  ' ms')
  cat ("\n   Number of Samples:     ", x$numSamples, '\n')
  cat ("\n")
}
