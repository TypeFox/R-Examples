##' collapse/bind several hyperSpec objects into one object
##'
##' The spectra from all objects will be put into one object.
##' The resulting object has all wavelengths that occur in the input objects.
##' Data points corresponding to wavelengths not in the original spectrum will be set to NA.
##' Extra data is combined in the same manner.
##' 
##' If the objects are named, the names will be preserved in extra data column \code{$.name}.
##' 
##' @author C. Beleites
##' @title Collapse hyperSpec objects
##' @export
##' @param ... hyperSpec objects to be collapsed into one object. Instead of giving several
##' arguments, a list with all objects to be collapsed may be given.
##' @aliases collapse collapse.hyperSpec
##' @seealso \code{\link[base]{merge}} to merge hyperSpec objects that share wavelengths but contain 
##'   different spectra,  \code{\link[base]{rbind}}, and  \code{\link[plyr]{rbind.fill}} for 
##' @return a hyperSpec object
##' @keywords manip
##' @examples
##' barbiturates [1:3]
##' barb <- collapse (barbiturates [1:3])
##' barb
##' 
##' a <- barbiturates [[1]]
##' b <- barbiturates [[2]]
##' c <- barbiturates [[3]]
##' 
##' a
##' b
##' c
##' collapse (a, b, c)
##' 

collapse <- function (...){
  dots <- list (...)

  ## accept also a list of hyperSpec objects
  if (length (dots) == 1 && is.list (dots [[1]]))
    dots <- dots [[1]]

  ## check the arguments
  lapply (dots, chk.hy)
  lapply (dots, validObject)

  ## names cause problems with unlisting labels.
  ## preserve them in column .name
  if (! is.null (names (dots))){
    dots <- mapply (function (object, name) {object$.name <- name; object}, dots, names (dots))
    names (dots) <- NULL
  }  
  
  ## prepare new labels
  labels <- unlist (lapply (dots, slot, "label"))
  labels <- labels [unique (names (labels))]
  
  ## merge data & spectra matrices
  dots <- lapply (dots, .wl2cln)
  dots <- rbind.fill (lapply (dots, slot, "data"))
  wl <- as.numeric (colnames (dots$spc))

  ## make a new hyperSpec object
  x <- new ("hyperSpec", wavelength = wl, data = dots, labels = labels)
  
  x
}

.wl2cln <- function (x){
  colnames (x@data$spc) <- formatC (x@wavelength, digits = 17)
  x
}

.test (collapse) <- function () {
  ## collapse messed up labels if a named list is collapsed
  tmp <- collapse (a = flu, b = flu)
  flu.labels <- lapply (flu@label, as.expression)
  checkEquals (labels (tmp) [names (flu.labels)], flu.labels)

  ## named lists should return .name column
  checkEquals (tmp$.name, rep (c ("a", "b"), each = nrow (flu)))
  
  ## no difference whether list or single arguments are given
  tmp2 <- list (a = flu, b = flu)
  tmp2 <- collapse (a = flu, b = flu)
  checkEquals (tmp, tmp2, 
               check.attributes = TRUE, check.names = TRUE, check.column.order = FALSE, check.label = TRUE)
}

## FIXME: 