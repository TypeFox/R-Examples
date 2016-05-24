# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

phasor = function (num, scaled = TRUE, add = FALSE, circle = FALSE, xlim,ylim, ...){
  if (!is.complex(num)) stop ('Input must be complex valued.')

  num = matrix(cbind(Re(num), Im(num)), length(num), 2)

  if (!scaled & missing(xlim)) xlim = max(abs(num[,1]))*c(-1,1) 
  if (!scaled & missing(ylim)) ylim = max(abs(num[,2]))*c(-1,1) 
  for (i in 1:nrow(num)){
    tmpnum = num[i,]
    if (scaled) tmpnum = tmpnum / sqrt(tmpnum[1]^2+tmpnum[2]^2)
    
	if (add & (dev.cur() == 1)) add = FALSE
    if (add | i > 1) arrows (0,0, tmpnum[1], tmpnum[2], ...)
    if (!add & i == 1){
      if (scaled) plot (0,0,type='n',xlim = c(-1.2,1.2), ylim = c(-1.2,1.2),xlab='Real Part', ylab='Imaginary Part')
      if (!scaled) plot (0,0,type='n',xlab='Real Part', ylab='Imaginary Part',xlim=xlim,ylim=ylim, ...)
      arrows (0,0, tmpnum[1], tmpnum[2], ...)
      if (circle) abline (h = 0, v = 0, lty = 'dotted')
      if (circle & scaled) sdellipse (matrix (c(1,0,0,1),2,2), means = c(0,0), stdev = 1)
    }
  }
}


