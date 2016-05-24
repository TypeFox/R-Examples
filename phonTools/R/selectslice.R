# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

selectslice = function (specobject, n = 1, plot = TRUE,...){
  if (!(class(specobject) == 'spectrogram')) stop ('Spectrogram object must be provided.')
  if (n < 1 | n %% 1 > 0) stop ('Positive integer n only.')

  if (plot) plot (specobject, ...)
  spect = specobject$spectrogram
  freqs1 = as.numeric(colnames(spect))
  times1 = as.numeric(rownames(spect))

  times = rep(times1, length(freqs1))
  freqs = rep(freqs1, each = length(times1))

  tmp = identify(times, freqs, "", n = n)
  time = sort(times[tmp])

  slices = NULL
  for (i in 1:n) slices = cbind (slices, spect[rownames(spect) == time[i],])
  colnames (slices) = time
  slices
}
