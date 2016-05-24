

splPrint <- function(x, ...){
  print.default(x, max = 10)
}

#' Printing spl object
#'
#' Print spl objects
#'
#'@param x an \code{spl} object.
#'@param ... additional arguments. Not used.
#'
#'@rdname printspl
#'@method print splpoints
#'@export
print.splpoints <- splPrint

#'@rdname printspl
#'@method print splicons
#'@export
print.splicons <- splPrint

#'@rdname printspl
#'@method print spllines
#'@export
print.spllines <- splPrint

#'@rdname printspl
#'@method print splpolygons
#'@export
print.splpolygons <- splPrint

#'@rdname printspl
#'@method print splgrid
#'@export
print.splgrid <- splPrint




#' Summary of map elements
#'
#' Get a summary of a map element.
#'
#'@param object a map layer, \code{basemap}, \code{spl} or \code{ui} object.
#'@param ... additional arguments. Not used.
#'
#'@rdname summarymaplayer
#'@method summary basemap
#'@export
summary.basemap <- function(object, ...){
  cat(object$name, ":\nrleafmap tile layer\nSource : ", object$URL)
}

#'@rdname summarymaplayer
#'@method summary splpoints
#'@export
summary.splpoints <- function(object, ...){
  cat(object$name, ":\nrleafmap vector layer of", length(object$coords[,1]), "points symbolized by circles \n")
}

#'@rdname summarymaplayer
#'@method summary splicons
#'@export
summary.splicons <- function(object, ...){
  cat(object$name, ":\nrleafmap vector layer of", length(object$coords[,1]), "points symbolized by icons \n")
}

#'@rdname summarymaplayer
#'@method summary spllines
#'@export
summary.spllines <- function(object, ...){
  cat(object$name, ":\nrleafmap vector layer of", length(object$coords), "lines \n")
}

#'@rdname summarymaplayer
#'@method summary splpolygons
#'@export
summary.splpolygons <- function(object, ...){
  cat(object$name, ":\nrleafmap vector layer of", length(object$coords), "polygons \n")
}

#'@rdname summarymaplayer
#'@method summary splgrid
#'@export
summary.splgrid <- function(object, ...){
  cat(object$name, ":\nrleafmap raster layer from a sp SpatialGridDataFrame\n------------\n")
  print(object$x)
  cat("------------")
}

#'@rdname summarymaplayer
#'@method summary ui
#'@export
summary.ui <- function(object, ...){
  cat("User interface configuration for rleafmap maps\n------------\n")
  cat("Zoom position :                ",object$zoom, "\n")
  cat("Layers control position :      ",object$layers, "\n")
  cat("Attributions position :        ",object$attrib, "\n")
  cat("Additionnal attribution text : ",object$attrib.text, "\n")
  cat("------------")
}