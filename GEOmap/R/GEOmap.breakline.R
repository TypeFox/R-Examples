GEOmap.breakline<-function(Z, ww)
{
  ########   ww is a vector of locations where there are
  #######   NA
  
  nww = length(ww)

  if(length(Z$x)<0)
    {
      x = Z[[1]]
      y = Z[[2]]
    }
  else
    {
           x = Z$x
           y = Z$y
    }

  newx = list()
  newy = list()

  if(nww<1) return(list(newx=x, newy=y) )

  ## ww = c(0, ww, length(Z$x)+1)
####  build up a set of strokes separated by NA's
  
  j = 1
  for(i in seq(from=1, to=nww, by=1)  )
    {

      i2 = (ww[i])

      
      newx[[i]] = x[j:i2]
      newy[[i]] =  y[j:i2]
      j = ww[i]+1
    }

  i2 = length(x)
  j=ww[i]+1
  i = i+1
  newx[[i]] = x[j:i2]
  newy[[i]] = y[j:i2]


  return(list(newx=newx, newy=newy))

}
