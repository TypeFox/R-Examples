sanity.check.fevdplot <-
function(fevddraws=fevddraws, type=type, labels=labels,save=save, bands=bands, grid=grid, bw=bw, table=table){
    #fevddraws
    Idim <- is.null(dim(fevddraws))

    if(Idim==TRUE){

      message(" ")
      stop(" fevddraws must be of dimensions (draws x steps x nvar) .\n", call. = FALSE)
    }

    Idiml <- length(dim(fevddraws))

    if(Idiml!=3){

      message(" ")
      stop(" fevddraws must be of dimensions (draws x steps x nvar) .\n", call. = FALSE)
    }

    Ina <- any(is.na(fevddraws)==TRUE)
    Inan <- any(is.nan(fevddraws)==TRUE)
    Inum <-  is.numeric(fevddraws)

    if(Ina==TRUE | Inan==TRUE){

      message(" ")
      stop(" fevddraws must not contain missing values.\n", call. = FALSE)
    }

    #type
    if(type!="median" & type!="mean"){

      message(" ")
      stop("FEVD type must be mean or median.\n", call. = FALSE)
    }
    #labels
    lbll <- length(labels)

    if(lbll < dim(fevddraws)[3]){

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
    return()
  }
