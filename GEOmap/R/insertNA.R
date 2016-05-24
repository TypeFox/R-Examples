insertNA<-function(y, ind)
  {
    wdiy = c(0, ind, length(y))
    j = 1
    dex = seq(from=(wdiy[j]+1), to=(wdiy[j+1]), by=1)
    
    Z = y[dex]
    for(j in 2:(length(wdiy)-1))
      {
        if((wdiy[j+1])<(wdiy[j]+1)) next
        dex = seq(from=(wdiy[j]+1), to=(wdiy[j+1]), by=1)
#### print(dex)
#### print(y[dex])
        Z = c(Z, NA, y[dex])
      }
    
    
    return(Z)
  }
