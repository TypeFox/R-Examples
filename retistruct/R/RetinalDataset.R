##' Construct an RetinalDataset that contains information specific
##' to the dataset in question. 
##'
##' @title RetinalDataset constructor
##' @param d A \code{dataset} object
##' @return An \code{retinalDataset} object. This contains all the
##' information from \code{d} plus:
##' \item{\code{DVflip}}{\code{TRUE} if the raw data is flipped in
##' the dorsoventral direction} 
##' \item{\code{side}}{The side of the eye ("Left" or "Right")}
##' @author David Sterratt
RetinalDataset <- function(d) {
  a <- d
  class(a) <- addClass("retinalDataset", a)
  a$DVflip <- FALSE
  a$side <- "Right"
  return(a)
}

##' Plot an retinal dataset. This basically is equivalent to
##' plotting a \code{dataset}, but may perform some transformations to
##' the date or plotting parameters. At present, if \code{DVflip} is
##' \code{TRUE}, it flips the \emph{y}-axis.
##'
##' @title Flat plot of retinal dataset
##' @param x \code{retinalDataset} object
##' @param axt whether to plot axes
##' @param ylim y-limits
##' @param ... Other plotting parameters
##' @method flatplot retinalDataset
##' @author David Sterratt
##' @export
flatplot.retinalDataset <- function(x, axt="n", ylim=NULL, ...) {
  if (x$DVflip) {
    if (is.null(ylim)) {
      ylim <- c(max(x$P[,2]), min(x$P[,2]))
    } else {
      ylim <- sort(ylim, TRUE)
    }
  }
  NextMethod(ylim=ylim)
}

##' Name a landmark in a \code{\link{RetinalDataset}}. This does the
##' same as the standard \code{\link{nameLandmark}}, but in addition,
##' if there exists a landmark named "OD", this creates a set of
##' points labelled "OD".
##'
##' @title Name a landmark in a RetinalDataset
##' @param d \code{\link{RetinalDataset}} object
##' @param i index of landmark to name
##' @param name name to give landmark
##' @return New \code{\link{RetinalDataset}} object in which landmark
##' is named
##' @author David Sterratt
##' @method nameLandmark retinalDataset
##' @export
nameLandmark.retinalDataset <- function(d, i, name) {
  d <- NextMethod()
  if (!is.na(getLandmarkID(d, "OD"))) {
    d$Ds[["OD"]] <- with(d, matrix(colMeans(Ss[[getLandmarkID(d, "OD")]]), 1, 2))
  }
  return(d)
}
