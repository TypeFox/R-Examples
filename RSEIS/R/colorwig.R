colorwig<-function(x1, y1, COL=rainbow(100))
{
  #### draw a time series with colors spanning te palette
  ###  this code uses segments to draw each line so it is
  ###   slow
  if(missing(COL)) COL=rainbow(100)
  nlen = length(x1)


  ncol = length(COL)

  KR = nlen/(ncol-1)
  
  KX = floor(seq(from=1, length=nlen)/(KR))+1
  
  cols=COL[KX]



  plot(x1, y1, type='n', xlab="time, s", ylab="Pa" )


  abline(h=0)
  segments(x1[1:(nlen-1)] , y1[1:(nlen-1)], x1[2:nlen], y1[2:nlen], col=cols)
  invisible(cols)	
}

