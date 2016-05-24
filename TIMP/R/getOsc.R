"getOsc" <-
  function (model) 
  {
    ## type is "harmonic"
    numosccol <- 0 
    if(tolower(model@oscspec$type) == "harmonic") {
      numosccol <- floor(length(model@oscpar)/3) 
      if (length(model@oscspec$start)>0 && length(model@oscpar)==0) {
        model@oscpar <- model@oscspec$start
      }
      if(length(model@nl)==0) {
        model@ncolc <- 
          array(model@ncomp + numosccol, 1) 
      } else {
        model@ncolc <- 
          array(model@ncomp + numosccol, model@nl) 
      }
    }    
    if(numosccol > 0) {
      model@cohcol <- (model@ncolc[1]-numosccol+1):model@ncolc[1]
    }
    
    model@ncomp <- max(model@ncolc)
    model
    
  }

