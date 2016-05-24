replaceWPX<-function(WPX, onepx , ind=1)
  {
    ####  delete members from  WPX list

    dpx = data.frame(WPX, stringsAsFactors = FALSE)
    opx = data.frame(onepx, stringsAsFactors = FALSE)
    
    dpx[ind,] = opx[1,]

    WPX = as.list(dpx)

    return(WPX)
  }
