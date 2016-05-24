# Copyright (c) 2015 Santiago Barreda
# All rights reserved.
spectralslice = function (sound, padding = length(sound) * 2, fs = 1, show = TRUE, 
    add = FALSE, window = "kaiser", windowparameter = 3, zeromax = TRUE, 
    preemphasisf = 0, type, line = FALSE, removeDC = TRUE, ...){
	
   if (class(sound) == "ts") fs = frequency (sound)
   if (class(sound) == "sound") {
        fs = sound$fs
        sound = sound$sound
    }
	
    if (preemphasisf > 0) 
        sound = preemphasis(sound, preemphasisf, fs)
    n = length(sound)
    if (removeDC) 
        sound = sound - mean(sound)
    sound = sound * windowfunc(n, window, windowparameter)
    N = n + padding
    if (fs > 1) 
        hz = seq(0, fs, length.out = N + 1)
    if (fs == 1) 
        hz = seq(0, 1, length.out = N + 1)
    hz = hz[hz <= fs/2]
    sound = c(sound, rep(0, padding))
    power = abs(fft(sound))
    power = power[1:length(hz)]/(n/2)
    power = log(power, 10) * 20
    power[which(power == min(power))] = sort(power)[2]
    if (zeromax == TRUE) 
        power = power - max(power)
    if (missing(type)) 
        type = "l"
    if (fs > 1) 
        xlab = "Frequency (Hz)"
    if (fs == 1) 
        xlab = "Frequency / Sampling Freq."
    if (add == FALSE & show == TRUE & line == FALSE) 
        plot(hz, power, ylab = "Power (dB)", xlab = xlab, type = type, 
            xaxs = "i", ...)
    if (add == TRUE & show == TRUE & line == FALSE) 
        lines(hz, power, type = type, ...)
    if (line == TRUE) {
        plot(hz, power, ylab = "Power (dB)", xlab = xlab, type = "p", 
            pch = 16, xaxs = "i", ...)
        segments(hz, rep(-5000, length(hz)), hz, power)
    }
    dB = power
    invisible(cbind(hz, dB))
}


slice = function (sound, padding = length(sound) * 2, fs = 1, show = TRUE, 
    add = FALSE, window = "kaiser", windowparameter = 3, zeromax = TRUE, 
    preemphasisf = 0, type, line = FALSE, removeDC = TRUE, ...){
  cl = match.call()
  args = sapply (2:length(cl), function(x) cl[[x]])
  names(args) = names(cl)[-1]
  do.call (spectralslice, args)
}

