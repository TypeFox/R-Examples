##' Create an object that is specific to retinal datasets. This
##' contains methods that return datapoint and landmark coordinates
##' that have been transformed according to the values of
##' \code{DVflip} and \code{side}.
##'
##' @title RetinalReconstructedDataset constructor
##' @param r Object that inherits both \code{reconstructedDataset} and
##' \code{dataset}.
##' @param report Function used to report progress.
##' @return \code{\link{RetinalReconstructedDataset}} object. This does not
##' contain any extra fields, but there are extra mthods dthat apply
##' to it.
##' @author David Sterratt
##' @export
RetinalReconstructedDataset <- function(r, report=message) {
  if (!(inherits(r, "reconstructedDataset") &
        inherits(r, "retinalDataset"))) {
    stop("Argument needs to inherit reconstructedDataset and retinalDataset")
  }
  class(r) <- addClass("retinalReconstructedDataset", r)
  report("Computing KDE and contours of points")
  r$KDE <- getKDE(r)
  report("Computing KR and contours of counts")
  r$KR <-  getKR(r)
  return(r)
}

##' Get spherical coordinates of datapoints, transformed according to
##' the values of \code{DVflip} and \code{side}.
##'
##' @title Get transformed spherical coordinates of datapoints
##' @param r \code{\link{RetinalReconstructedDataset}} object.
##' @return \code{Dss}
##' @method getDss retinalReconstructedDataset
##' @author David Sterratt
getDss.retinalReconstructedDataset <- function(r) {
  Dss <- NextMethod()
  if (length(Dss) > 0) {
    if (r$DVflip) {
      for (i in 1:length(Dss)) {
        Dss[[i]][,"lambda"] <- -Dss[[i]][,"lambda"]
      }
    }
  }
  return(Dss)
}


##' @title Get grouped variable with locations in spherical coordinates.
##' @param r \code{\link{ReconstructedDataset}} or \code{\link{RetinalReconstructedDataset}} object.
##' @return \code{Gss}
##' @method getGss retinalReconstructedDataset
##' @author David Sterratt
##' @export
getGss.retinalReconstructedDataset <- function(r) {
  Gss <- NextMethod()
  if (length(Gss) > 0) {
    if (r$DVflip) {
      for (i in 1:length(Gss)) {
        Gss[[i]][,"lambda"] <- -Gss[[i]][,"lambda"]
      }
    }
  }
  return(Gss)
}

##' Get Karcher mean of datapoints in spherical coordinates,
##' transformed according to the values of \code{DVflip} and
##' \code{side}.
##'
##' @title Get transformed spherical coordinates of Karcher mean of
##' datapoints
##' @param r \code{\link{RetinalReconstructedDataset}} object.
##' @return \code{Dss.mean}
##' @method getDssMean retinalReconstructedDataset
##' @author David Sterratt
getDssMean.retinalReconstructedDataset <- function(r) {
  Dss.mean <- NextMethod()
  Dss.mean[["OD"]] <- NULL
  if (length(Dss.mean) > 0) {
    if (r$DVflip) {
      for (i in 1:length(Dss.mean)) {
        Dss.mean[[i]][,"lambda"] <- -Dss.mean[[i]][,"lambda"]
      }
    }
  }
  return(Dss.mean)
}

##' Get spherical coordinates of landmarks, transformed according to
##' the values of \code{DVflip} and \code{side}.
##'
##' @title Get transformed spherical coordinates of datapoints
##' @param r \code{\link{RetinalReconstructedDataset}} object.
##' @return \code{Dss}
##' @method getSss retinalReconstructedDataset
##' @author David Sterratt
getSss.retinalReconstructedDataset <- function(r) {
  Sss <- NextMethod()
  if (length(Sss) > 0) {
    if (r$DVflip) {
      for (i in 1:length(Sss)) {
        Sss[[i]][,"lambda"] <- -Sss[[i]][,"lambda"]
      }
    }
  }
  return(Sss)
}

##' @export
getSssMean.retinalReconstructedDataset <- function(r) {
  Sss.mean <- NextMethod()
  if (length(Sss.mean) > 0) {
    if (r$DVflip) {
      for (i in 1:length(Sss.mean)) {
        Sss.mean[[i]][,"lambda"] <- -Sss.mean[[i]][,"lambda"]
      }
    }
  }
  return(Sss.mean)
}

##' @method projection retinalReconstructedDataset
##' @export
projection.retinalReconstructedDataset <- function(r, ...) {
  ## This will call projection.reconstructedDataset()
  NextMethod()
}


