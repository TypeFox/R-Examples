# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

plot.spectrogram = function (x, y, ylim,xlim, quality = FALSE, ...){
  if (x$colors[1] == TRUE)
    zcolors = colorRampPalette (c('dark blue','blue','cyan','light green','yellow',
                                  'orange','red', 'brown'))
  if (x$colors[1] == FALSE) zcolors = colorRampPalette (c('white','black'))
  if (length(x$colors) > 1) zcolors = colorRampPalette (x$colors)
  
  zrange = c(-x$dynamicrange,0)
  nlevels = abs (zrange[1] - zrange[2]) * 1.2
  
  levels = pretty(zrange, nlevels)
  zcolors = zcolors(length(levels) - 1);
  times = as.numeric(rownames (x$spectrogram))
  hz = as.numeric(colnames (x$spectrogram))

  x$spectrogram[which(x$spectrogram < (-1 * x$dynamicrange))] = -1 * x$dynamicrange
  
  if (missing (ylim)) ylim = range (0, x$maxfreq)
  if (missing (xlim)) xlim = range (times)
 
  if (quality){  
    plot.new()
    plot.window(ylim = ylim, xlim = xlim, xaxs = 'i', yaxs = 'i', ...)
    .filled.contour(as.double(times), as.double(hz), x$spectrogram, as.double(levels), col = zcolors)
    Axis(times, side = 1)
    Axis(hz, side = 2)
    title(xlab = "Time (ms)", ylab = "Frequency (Hz)")
  }
  if (!quality){
    image (as.double(times), as.double(hz), x$spectrogram, useRaster = FALSE, col = zcolors, xlab = 'Time (ms)', ylab = "Frequency (Hz)", ylim = ylim, xlim = xlim,...)
 }
 box()
}

