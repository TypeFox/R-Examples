sanity.check.reject <-
function(Y=Y, nlags=nlags, draws=draws, subdraws=subdraws, nkeep=nkeep, KMIN=KMIN, KMAX=KMAX, constrained=constrained, constant=constant, steps=steps){

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


      subdrawsint <- subdraws%%1==0

      if(subdrawsint!=TRUE){

        message(" ")
        stop(" Number of subdraws must be integer.\n", call. = FALSE)
    }

      if(subdraws<=0){

        message(" ")
        stop(" Number of subdraws must geater than zero.\n", call. = FALSE)
    }

      nkeepint <- nkeep%%1==0

      if(nkeepint!=TRUE){

        message(" ")
        stop(" Number of nkeep must be integer.\n", call. = FALSE)
    }

      if(nkeep<=0){

        message(" ")
        stop(" Number of nkeep must geater than zero.\n", call. = FALSE)
    }


      KMINint <- KMIN%%1==0

      if(KMINint!=TRUE){

        message(" ")
        stop("KMIN must be integer.\n", call. = FALSE)
    }

      if(KMIN<=0 | KMIN>KMAX){

        message(" ")
        stop("KMIN must be greater than zero and smaller than KMAX.\n", call. = FALSE)
      }


      KMAXint <- KMAX%%1==0

      if(KMAXint!=TRUE){

        message(" ")
        stop("KMAX must be integer.\n", call. = FALSE)
    }

      if(KMIN<=0 | KMIN>KMAX){

        message(" ")
        stop("KMAX must be greater than zero and greater than KMIN.\n", call. = FALSE)
      }


 cnsl <- length(constrained)

    if(cnsl<=0){

      message(" ")
        stop("Number of constraints must at least 1.\n", call. = FALSE)
    }

        if(cnsl>ncol(Y)){

      message(" ")
        stop("Number of constraints cannot be greater than number of variables.\n", call. = FALSE)
    }

   if(max(abs(constrained))>ncol(Y) | min(abs(constrained))==0){

      message(" ")
        stop("Constraints must be between 1 and the number of variables.\n", call. = FALSE)
    }

 cnscns <- all(constrained%%1==0)

    if(cnscns==FALSE){

      message(" ")
        stop("All constraints must be integers.\n", call. = FALSE)
 }

cnsdup <- anyDuplicated(abs(constrained))
    if(cnsdup>0){

      message(" ")
        stop("Cannot provide multiple constraints for the same variable.\n", call. = FALSE)
 }

 return()

}
