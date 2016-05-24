sanity.check.target <-
function(Y=Y, nlags=nlags, irfdraws=irfdraws, constant=constant, type=type, labels=labels, target= target, save=save, legend=legend, bands=bands, grid=grid, bw=bw, maxit=maxit){


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


  Idim <- is.null(dim(irfdraws))

  if(Idim==TRUE){

    message(" ")
    stop(" Irfdraws must be of dimensions (draws x steps x nvar) .\n", call. = FALSE)
  }

  Idiml <- length(dim(irfdraws))

  if(Idiml!=3){

    message(" ")
    stop(" Irfdraws must be of dimensions (draws x steps x nvar) .\n", call. = FALSE)
  }

  Ina <- any(is.na(irfdraws)==TRUE)
  Inan <- any(is.nan(irfdraws)==TRUE)
  Inum <-  is.numeric(irfdraws)

  if(Ina==TRUE | Inan==TRUE){

    message(" ")
    stop(" Irfdraws must not contain missing values.\n", call. = FALSE)
  }

  #type
  if(type!="median" & type!="mean"){

    message(" ")
    stop("IRF type must be mean or median.\n", call. = FALSE)
  }
  #labels
  lbll <- length(labels)

  if(lbll < dim(irfdraws)[3]){

    message(" ")
    stop("Number of labels must be equal number of variables in the model.\n", call. = FALSE)
  }
  bndtest <- is.null(bands)
  if(bndtest!=TRUE){
    bndsl <- length(bands)

    if(bndsl !=2){

      message(" ")
      stop("Error bands must contain only two values c(lower, upper) or 'NULL'.\n", call. = FALSE)
    }


    if(max(bands)>=1 | min(bands)<=0){

      message(" ")
      stop("Error bands must be between 0 and 1.\n", call. = FALSE)
    }
  }




mxint <- maxit%%1==0

if(mxint!=TRUE){

  message(" ")
  stop(" Number of maxit must be integer.\n", call. = FALSE)
}

if(maxit<=0){

  message(" ")
  stop(" Number of maxit must geater than zero.\n", call. = FALSE)
}
  return()
}
