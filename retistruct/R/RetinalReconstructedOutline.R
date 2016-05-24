##' Create an object that is specific to retinal datasets. This
##' contains methods that return datapoint and landmark coordinates
##' that have been transformed according to the values of
##' \code{DVflip} and \code{side}.
##'
##' @title RetinalReconstructedOutline constructor
##' @param r Object that inherits \code{ReconstructedOutline}
##' @param report Function used to report progress.
##' @return \code{RetinalReconstructedOutline} object. This does not
##' contain any extra fields, but there are extra mthods dthat apply
##' to it.
##' @author David Sterratt
##' @export
RetinalReconstructedOutline <- function(r, report=message) {
  if (!(inherits(r, "reconstructedOutline"))) {
    stop("Argument needs to inherit reconstructedOutline")
  }
  class(r) <- addClass("retinalReconstructedOutline", r)
  return(r)
}

##' @method getIms retinalReconstructedOutline
##' @export
getIms.retinalReconstructedOutline <- function(r) {
  ims <- NextMethod()
  if (r$DVflip) {
    if (!is.null(ims)) {
      ims[,"lambda"] <- -ims[,"lambda"]
    }
  }
  return(ims)
}

##' @export
getTss.retinalReconstructedOutline <- function(r) {
  Tss <- NextMethod()
  if (r$DVflip) {
    for (i in 1:length(Tss)) {
      Tss[[i]][,"lambda"] <- -Tss[[i]][,"lambda"]
    }
  }
  return(Tss)
}

##' @method projection retinalReconstructedOutline
##' @export
projection.retinalReconstructedOutline <-
  function(r,
           projection=azimuthal.equalarea,
           ...) {
  philim <- c(-90, 90)
  colatitude <- FALSE
  pole <- TRUE
  if (!(identical(projection, sinusoidal) |
        identical(projection, orthographic))) {
    philim <- c(-90, r$phi0*180/pi)
    colatitude <- TRUE
    pole <- FALSE
  }
  if (r$side=="Right") {
    labels=c("N", "D", "T", "V")
  } else {
    labels=c("T", "D", "N", "V")
  }
  NextMethod(philim=philim,
             labels=labels,
             colatitude=TRUE)
}
