sanity.check.irfplot <-
function(irfdraws=irfdraws,type=type, labels=labels,save=save, bands=bands, grid=grid, bw=bw){
#irfdraws
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
    return()
}
