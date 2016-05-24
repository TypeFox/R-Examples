`lagplot` <-
function(y1,dt, lag, PLOT=FALSE)
{
  if(missing(PLOT)) { PLOT=TRUE }
  ##  y1 = input time series
  ##  dt deltaT for time series
  ##  lag in time for y1

  N = length(y1)
  ex = seq(from=0, length=N, by = 0.008)

  plot(ex, y1, type='l')
  nlag = round(lag/dt)

  if(nlag>0)
    {
      y3 =  c(rep(NA, length=abs(nlag)), y1[1:(N-nlag)])

    }
  else
    {
      y3 =  c( y1[(nlag):N], rep(NA, length(abs(nlag))))
    }

  if(PLOT)
    {
      plot(ex,y3, type='l')
    }
  invisible(y3)

}

