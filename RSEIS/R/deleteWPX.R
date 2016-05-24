deleteWPX<-function(WPX, ind=1)
  {
    ####  delete members from  WPX list

    dpx = data.frame(WPX, stringsAsFactors = FALSE)
    
    dpx = dpx[-ind,]

    WPX = as.list(dpx)

    return(WPX)
  }
