NadarayaWatsonkernel <- function(x, y, h, gridpoint)
{
  n = length(y)
  mh = vector(,length(gridpoint))
  for(j in 1:length(gridpoint))
  {
      suma = sumb = vector(,n)
      for(i in 1:n)
      {
          suma[i] = ker((gridpoint[j] - x[i])/h) * y[i]
          sumb[i] = ker((gridpoint[j] - x[i])/h)
      }
      mh[j] = sum(suma)/sum(sumb)
  }
  return(list(gridpoint = gridpoint, mh = mh))
}

