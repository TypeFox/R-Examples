GEOmap.breakpoly<-function(Z, ww)
{
  ########   ww is a list of locations where there should be breaks
  #######  
  nww = length(ww)
  newx = list()
  newy = list()
  if(nww<1) return(list(newx=Z$x, newy=Z$y) )
  ##
####  build up a set of strokes separated by NA's
  ##  connect the last stroke to the first

  i2 = length(Z$x)
  j=ww[nww]+1
  i = 1
  newx[[1]] = c(   Z$x[j:i2], Z$x[ 1:(ww[1]) ])
  newy[[1]] = c(   Z$y[j:i2], Z$y[ 1:(ww[1]) ])


  
  j = ww[1]+1
  for(i in seq(from=2, to=nww, by=1)  )
    {
      i2 = (ww[i])
      newx[[i]] = Z$x[j:i2]
      newy[[i]] =  Z$y[j:i2]
      j = ww[i]+1
    }


  return(list(newx=newx, newy=newy))

}
