
#' Defined generics
#'
#' Generics defined by the archetypes package.
#'
#' @param object An object
#' @param ... Futher arguments
#' @rdname archetypes-generics
#'
#' @export
rss <- function(object, ...) {
  UseMethod('rss')
}



#' Return number of parameters
#'
#' @rdname archetypes-generics
#'
#' @export
nparameters <- function(object, ...) {
  UseMethod('nparameters')
}




#' Return best model
#'
#' @rdname archetypes-generics
#'
#' @export
bestModel <- function(object, ...) {
  UseMethod('bestModel')
}



#' Panorama
#'
#' @rdname archetypes-generics
#'
#' @export
panorama <- function(object, ...) {
  UseMethod('panorama')
}



#' Parallel coordinates plot
#'
#' @param x An object.
#' @rdname archetypes-generics 
#'
#' @export
pcplot <- function(x, ...) {
  UseMethod('pcplot')
}

