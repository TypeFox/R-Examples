##' This adds a new class to the class vector. If the class is dataset
##' type it is prepended at the start of the class vector. If it is an
##' outline type, it is put after all the dataset classes, but before
##' the other outline ones. This is needed for the plotting functions
##' to work properly.
##'
##' @title Add new class to class vector of object
##' @param newclass New class to add
##' @param obj Object to which to add it
##' @return New class vector
##' @author David Sterratt
addClass <- function(newclass, obj) {
  if (newclass %in% class(obj)) {
    return(class(obj))
  }
  cl <- unique(c(newclass, class(obj)))
  cld <- cl[grep("[Dd]ataset", cl)]
  clo <- cl[setdiff(1:length(cl), grep("[Dd]ataset", cl))]
  return(c(cld, clo))
}

##' Plot flat representation of object
##'
##' @title Flat plot of object
##' @param x \code{\link{Outline}}, \code{\link{Dataset}} \&c object
##' @param axt whether to plot axes
##' @param ylim y-limits
##' @param ... Other plotting parameters
##' @author David Sterratt
##' @export
flatplot <- function(x, axt="n", ylim=NULL, ...) {
  UseMethod("flatplot")
}

##' @export
flatplot.default <- function(x, axt="n", ylim=NULL, ...) {
}

##' Generic function for plotting projections of objects.
##'
##' @title Plot projection of an object
##' @param r Object such as a \code{\link{ReconstructedOutline}}
##' @param ... Other parameters; see
##' \code{\link{projection.reconstructedOutline}} and
##' \code{\link{projection.reconstructedDataset}}
##' @author David Sterratt
##' @export
projection <- function(r, ...) {
  UseMethod("projection")
}

##' @export
projection.default <- function(r, ...) {
  plot.new()
}

##' Draw a spherical plot.
##'
##' @title Spherical plot of reconstructed outline
##' @param r Object with class \code{\link{ReconstructedOutline}}
##' and/or \code{\link{ReconstructedDataset}}
##' @param ... Parameters depending on class of \code{r}
##' @author David Sterratt
##' @export
sphericalplot <- function(r, ...) {
  UseMethod("sphericalplot")
}

##' @export
sphericalplot.default <- function(r, ...) {
  rgl.clear()
  rgl.bg(color="white")
}


##' Plot the fractional change in length of mesh edges. The length of
##' each edge in the mesh in the reconstructed object is plotted
##' against each edge in the spherical object. The points are
##' colour-coded according to the amount of log strain each edge is
##' under.
##'
##' @title Plot the fractional change in length of mesh edges
##' @param r \code{\link{ReconstructedOutline}} object
##' @author David Sterratt
##' @export
lvsLplot <- function(r) {
  UseMethod("lvsLplot")
}

##' @export
lvsLplot.default <- function(r) {
}

##' Get spherical coordinates of datapoints.
##'
##' @title Get transformed spherical coordinates of datapoints
##' @param r \code{\link{ReconstructedDataset}} or \code{\link{RetinalReconstructedDataset}} object.
##' @return \code{Dss}
##' @author David Sterratt
##' @export
getDss <- function(r) {
  UseMethod("getDss")
}

##' @export
getDss.default <- function(r) {
  return(NULL)
}

##' Get Karcher mean of datapoints in spherical coordinates.
##'
##' @title Karcher mean of datapoints in spherical coordinates
##' @param r \code{\link{ReconstructedDataset}} or \code{\link{RetinalReconstructedDataset}} object.
##' @return \code{Dss.mean}
##' @author David Sterratt
##' @export
getDssMean <- function(r) {
  UseMethod("getDssMean")
}

##' @export
getDssMean.default <- function(r) {
  return(NULL)
}

##' @title Get grouped variable with locations in spherical coordinates.
##' @param r \code{\link{ReconstructedDataset}} or \code{\link{RetinalReconstructedDataset}} object.
##' @return \code{Gss}
##' @author David Sterratt
##' @export
getGss <- function(r) {
  UseMethod("getGss")
}

##' @export
getGss.default <- function(r) {
  return(NULL)
}

##' Get spherical coordinates of landmarks.
##'
##' @title Get transformed spherical coordinates of landmarks.
##' @param r \code{\link{ReconstructedDataset}} or \code{\link{RetinalReconstructedDataset}} object.
##' @return \code{Sss}
##' @author David Sterratt
##' @export
getSss <- function(r) {
  UseMethod("getSss")
}

##' @export
getSss.default <- function(r) {
  return(NULL)
}

##' Get Karcher mean of landmarks in spherical coordinates. If applied
##' to a \code{\link{RetinalReconstructedDataset}}, the coordinates
##' are flipped if \code{DVflip} is \code{TRUE}.
##'
##' @title Karcher mean of landmarks in spherical coordinates
##' @param r \code{\link{ReconstructedDataset}} or \code{\link{RetinalReconstructedDataset}} object.
##' @return List of spherical coordinates, with one named element per
##' landmark.
##' @author David Sterratt
##' @export
getSssMean <- function(r) {
  UseMethod("getSssMean")
}

##' @export
getSssMean.default <- function(r) {
  return(NULL)
}

##' Get spherical coordinates of tears. If applied to a
##' \code{\link{RetinalReconstructedOutline}}, the tears are flipped
##' according to whether the retina is DV flipped or not.
##'
##' @title Get spherical coordinates of tears.
##' @param r \code{\link{ReconstructedOutline}} or
##' \code{\link{RetinalReconstructedOutline}} object.
##' @return List of tears in spherical coordinates.
##' @author David Sterratt
##' @export
getTss <- function(r) {
  UseMethod("getTss")
}

##' @export
getTss.default <- function(r) {
  return(NULL)
}

##' @title Get coordinates of corners of pixels of image in spherical coordinates 
##' @param r \code{\link{ReconstructedOutline}} or
##' \code{\link{RetinalReconstructedOutline}} object
##' @return Coordinates of corners of pixels in spherical coordinates 
##' @author David Sterratt
##' @export
getIms <- function(r) {
  UseMethod("getIms")
}

##' @export
getIms.default <- function(r) {
  return(NULL)
}

##' @title Get IDs of groups of data within a dataset
##' @export
##' @param r Object
getIDs <- function(r) {
  UseMethod("getIDs")
}

##' @export
getIDs.default <- function(r) {
  return(NULL)
}

##' Name a landmark in a \code{\link{Dataset}}. The name of element
##' \code{i} of \code{Ss} is set to \code{name}, the name of any
##' element that bore the name is set to "" and all other elements are
##' unaltered.
##'
##' @title Name a landmark in a Dataset
##' @param d \code{\link{Dataset}} object
##' @param i index of landmark to name
##' @param name name to give landmark
##' @return New \code{\link{Dataset}} object in which landmark is
##' named
##' @author David Sterratt
##' @export
nameLandmark <- function(d, i, name) {
  UseMethod("nameLandmark")
}

##' @export
nameLandmark.default <- function(d, i, name) {
  return(d)
}


