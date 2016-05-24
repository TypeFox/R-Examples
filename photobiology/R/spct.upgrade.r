#' Upgrade one spectral object
#'
#' Update the spectral class names of objects to those used in photobiology (>=
#' 0.6.0) and add 'version' attribute as used in photobiology (>= 0.70).
#'
#' @param object generic.spct A single object to upgrade
#'
#' @note The object is modified by reference. The class names with ending
#'   ".spct" replaced by their new equivalents ending in "_spct".
#'
#' @return The modified object (invisibly).
#'
#' @export
#'
#' @family upgrade from earlier versions
#'
upgrade_spct <-
  function(object) {
    name <- substitute(object)
    version <- getSpctVersion(object)
    if (version < 1L) {
    class(object) <- gsub(".spct", "_spct", class(object), fixed = TRUE)
    }
    if (version <= 1L) {
      attr(object, "spct.version") <- 2L
    }
    check_spct(object)
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, object, parent.frame(), inherits = TRUE)
    }
    invisible(object)
  }

#' Upgrade one or more spectral objects
#'
#' Update the spectral class names of objects to those used in photobiology (>=
#' 0.6.0).
#'
#' @param obj.names char Names of objects to upgrade as a vector of character
#'   strings
#'
#' @note The objects are modified by reference. The class names with ending
#'   ".spct" are replaced by their new equivalents ending in "_spct".
#'   \code{object.names} can safely include names of any R object. Names of
#'   objects which do not belong to any the old \code{.spct} classes are
#'   ignored. This makes it possible to supply as argument the output from
#'   \code{ls}, the default, or its equivalent \code{objects}.
#'
#' @return The modified object (invisibly).
#'
#' @export
#'
#' @family upgrade from earlier versions
#'
upgrade_spectra <- function(obj.names = ls(parent.frame())) {
  for (obj.name in obj.names) {
    obj <- get(obj.name, inherits = TRUE)
    if (is.old_spct(obj)) {
      message("Upgrading: ", obj.name)
      upgrade_spct(obj)
      assign(obj.name, obj, inherits = TRUE)
    } else {
      message("Skipping: ", obj.name)
    }
  }
}

#' Query if an object has old class names
#'
#' Query if an object has old class names Query if an object has old class names
#' as used in photobiology (>= 0.6.0).
#'
#' @param object an R object
#'
#' @return logical
#'
#' @export
#'
#' @family upgrade from earlier versions
#'
is.old_spct <- function(object) inherits(object, "generic.spct")
