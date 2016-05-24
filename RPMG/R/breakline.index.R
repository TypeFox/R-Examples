breakline.index<-function(Z, ww)
{
  #######   Z is a vector of indices
  ########   ww is a vector of locations
  ########       where there should be breaks
  #######

  if(is.matrix(ww))
    {
      jays = ww[,1]
      endz = ww[,2]
    }
  if(is.vector(ww))
    {
      jays = c(1, ww+1) 
      endz = c(ww, length(Z)) 
    }

 ##### print(jays)
 ##### print(endz)
  
  nww = length(ww)

  if(nww<1) return(list(Z))
  x = Z
  newind  = list()
  
  nww = length(jays)
  
  for(i in seq(from=1, to=nww, by=1)  )
    {
      newind[[i]] = x[jays[i]:endz[i] ]
     
    }

  return(newind)
}
