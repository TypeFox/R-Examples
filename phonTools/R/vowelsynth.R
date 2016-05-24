# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

vowelsynth = function (ffs = c(270, 2200, 2800, 3400, 4400), fbw = 0.06, dur = 300, f0 = c(120,100), 
                       fs = 10000, verify = FALSE, returnsound = TRUE, noise1 = .001, noise2 = .01, power = NULL){
  if (dur < 50) stop ('Duration must be at least 50 ms.')
  if (length (f0) > 2) stop ('Only initial and final f0 may be specified.')
  
  T = 1/fs
  n = round(dur/1000/T)
  if (!is.null(power)) n = length(power)
  if (is.numeric (nrow(ffs))) n = nrow(ffs)

  if (length(f0) == 1) f0 = c(f0,f0)
  f0 = exp(seq (log(f0[1]), log(f0[2]), length.out = n))
 
  vsource = NULL
  spot = 1
  while (length (vsource) < n*5){
    tmp = f0[spot]
    cycle = round(fs/tmp)
    tmp = 2*seq (0,1, 1/(cycle*4)) - 3*seq (0,1, 1/(cycle*4))^2 
    tmp = c(rep (0, cycle), tmp)
    vsource = c(vsource, tmp)
    spot = spot + cycle
  }
  vsource = resample (vsource, fs, fs*5, synthfilter = TRUE)
  vsource = vsource + rnorm (length(vsource), sd = sd(vsource)*noise1)
  vsource = vsource[1:n] 
  
  vsource = jitter(vsource)
  vsource = preemphasis (vsource, coeff = .94)

  x = c(1, 10 / (1000/fs), 20 / (1000/fs), n-(30 / (1000/fs)), n)
  if (is.null(power))   power = interpolate (x, y = c(30,55, 60,55,30), increment = 1, type = 'linear')[1:n,2]
  power = 10^(power/20)
  
  power = jitter(power, factor = 0.01)
  vsource = vsource * power
  output = Ffilter (vsource, ffs = ffs, fs = fs, verify = FALSE, bwp = fbw)
  
  output = output * power
  output = output + rnorm (length(output), sd = sd(output)*noise2)
  output = output /(max(abs(output)) * 1.05)
  
  if (verify == TRUE) {
    par (mfrow = c(2,1), mar = c(4,4,1,1))
    plot ((1:n)*T*1000,output, ylab = 'Amplitude', xlab = 'Time (ms)', type = 'l', xaxs = 'i')
    abline (h = 0, lty = 'dotted') 
    spectrogram (output, fs = fs, dynamicrange = 60)
  }
  if (returnsound == TRUE) 
    output = makesound(output, "sound.wav", fs = fs)
  return(output)
}





