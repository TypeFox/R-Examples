 ridNA<-function(z, temp)
  {
    ######  replace NA with something else (numeric)
    z[is.na(z)] = temp
    return(z)
  }
  


