importance <- function(x, ...) UseMethod("importance")

importance.wsrf <- function(x, type=NULL, class=NULL, scale=TRUE, ...) {

  imp        <- x[[.IMPORTANCE_IDX]]
  impSD      <- x[[.IMPORTANCESD_IDX]]
  hasPermImp <- !is.null(impSD)
  igrIdx     <- ncol(imp)
  hasType    <- !is.null(type)
  hasClass   <- !is.null(class)
  
  ## For unmached arguments.

  if (!hasPermImp && ((hasType && type == 1) || (!hasType && hasClass)))
    stop("That measure has not been computed.")  # No Perm-based measures, but required.

  if (hasType) {
    if (type != 1 && type != 2) stop("Wrong type specified.")
    if (type == 2 && hasClass)  stop("No class-specific measure for that type.")
  }


  ## When arguments matched and result needs to be converted.

  if (hasType && type == 2)                      # Only MeanDecreaseIGR
    imp <- imp[, igrIdx, drop=FALSE]

  
  if (hasType && type == 1 && !hasClass) {       # Only MeanDecreaseAccuracy
    if (scale) imp <- imp[, -igrIdx] / impSD
    else       imp <- imp[, -igrIdx]
  }
  
  if (!hasType && !hasClass && scale)            # Scaled MeanDecreaseAccuracy and MeanDecreaseIGR
    imp <- cbind(imp[, -igrIdx] / impSD, imp[, igrIdx, drop=FALSE])


  if ((!hasType || type == 1) && hasClass) {     # Class-specific measures

    whichCol   <- match(class, colnames(imp)[-igrIdx])
    whichNA    <- which(is.na(whichCol))

    if (length(whichNA) != 0) {
      outmessage <- paste("Class", paste(class[whichNA], collapse=", "), "not found.")
      if (length(whichNA) != length(whichCol))
        warning(outmessage)  # Some classes matched.
      else
        stop(outmessage)     # No class matched.
    }

    if (length(whichNA)!=0) whichCol <- whichCol[-whichNA]

    if (scale) imp <- imp[, whichCol, drop=FALSE] / impSD[, whichCol, drop=FALSE]
    else       imp <- imp[, whichCol, drop=FALSE]
  }

  imp
}
