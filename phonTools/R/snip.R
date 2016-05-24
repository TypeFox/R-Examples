# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

snip = function (object, show = TRUE){
  if (class(object) != "spectrogram" & class(object) != "sound") 
    stop("Input must be a sound or spectrogram object")
  if (class(object) == "sound") {
    time = 1:length(object$sound)/object$fs * 1000
    y = object$sound
    touse = seq(1, length(time), 2)
    plot(time[touse], y[touse], xlab = "Time (ms)", ylab = "Amplitude", 
         type = "l", xaxs = "i")
    times = seq(min(time), max(time), length.out = 1000)
    amps = seq(min(y), max(y), length.out = 100)
    edges = identify(rep(times, length(amps)), rep(amps, 
                                                   each = length(times)), "", n = 2)
    edges = edges%%1000
    edges = sort(times[edges])
    T = 1/object$fs
    start = edges[1]/1000/T
    end = edges[2]/1000/T
    snipped = object$sound[start:end]
    newtime = T * 1000 * (1:length(snipped))
    newsound = makesound(snipped, object$filename, object$fs)
    if (show == TRUE) 
      plot(newtime[seq(1, length(snipped), 2)], snipped[seq(1,length(snipped), 2)], 
      xlab = "Time (ms)", ylab = "Amplitude", type = "l", xaxs = "i")
    output = newsound
  }
  if (class(object) == "spectrogram") {
    plot(object)
    times = as.numeric(rownames(object$spectrogram))
    times = seq(head(times, 1), tail(times, 1), length.out = 1000)
    freqs = as.numeric(colnames(object$spectrogram))
    freqs = seq(head(freqs, 1), tail(freqs, 1), length.out = 100)
    edges = identify(rep(times, length(freqs)), rep(freqs, 
                                                    each = length(times)), "", n = 2)
    edges = sort(edges%%1000)
    specttimes = as.numeric(rownames(object$spectrogram))
    object$spectrogram = object$spectrogram[specttimes > times[edges[1]] & specttimes < times[edges[2]], ]
    tmp = as.numeric(rownames(object$spectrogram))
    rownames(object$spectrogram) = tmp - min(tmp)
    if (show == TRUE)   plot(object)
    output = object
  }
  invisible(output)
}

