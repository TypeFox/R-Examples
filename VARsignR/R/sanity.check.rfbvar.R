sanity.check.rfbvar <-
function(Y=Y, nlags=nlags, draws=draws, constant=constant, steps=steps, shock=shock){

  Yts <-is.ts(Y)
  Ydf <- is.data.frame(Y)

if(Yts==FALSE & Ydf==FALSE){
  
  message(" ")
    stop(" Data has to be a data.frame() or ts() object.\n", call. = FALSE)
}

  if(ncol(Y)<2){
    
    message(" ")
    stop(" Need more than 1 variable.\n", call. = FALSE)
  }

    Yna <- any(is.na(Y)==TRUE)
    Ynan <- any(is.nan(Y)==TRUE)
    Ynum <-  is.numeric(Y)

    if(Yna==TRUE | Ynan==TRUE){
      
      message(" ")
        stop(" Data must not contain missing values.\n", call. = FALSE)
    }

    if(Ynum!=TRUE){
      
      message(" ")
        stop(" Data must not contain strings.\n", call. = FALSE)
    }

    nlagsint <- nlags%%1==0
    nlagsobs <- nrow(Y)-nlags

    if(nlagsint==FALSE){
      
      message(" ")
        stop("Number of lags must be integer.\n", call. = FALSE)
    }

      if(nlagsobs<=0){
        
        message(" ")
        stop(" Number of lags cannot be larger than the number of observations.\n", call. = FALSE)
    }

  if(nlags<1){
    
    message(" ")
    stop(" Need at least 1 lag.\n")
}
        if(nlags>nrow(Y)){
          
          message(" ")
    stop(" Number of lags have to be smaller than number of observations.\n", call. = FALSE)
  }

    drawsint <- draws%%1==0

      if(drawsint!=TRUE){
        
        message(" ")
        stop(" Number of draws must be integer.\n", call. = FALSE)
    }

      if(draws<=0){
        
        message(" ")
        stop(" Number of draws must geater than zero.\n", call. = FALSE)
    }

     stepsint <- steps%%1==0

      if(stepsint!=TRUE){
        
        message(" ")
        stop(" Number of steps must be integer.\n", call. = FALSE)
    }

      if(steps<=0){
        
        message(" ")
        stop(" Number of steps must geater than zero.\n", call. = FALSE)
    }

     shockint <- shock%%1==0

      if(shockint!=TRUE){
        
        message(" ")
        stop(" Shock of interest must be integer.\n", call. = FALSE)
    }

      if(shock<=0 | shock>ncol(Y)){
        
        message(" ")
        stop(" Shock of interest must be greater than zero and smaller or equal number variables.\n", call. = FALSE)
      }
 return()
}
