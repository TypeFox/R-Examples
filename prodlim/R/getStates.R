##' Extract the states of a multi-state model
##'
##' Applying this function to the fit of prodlim means to apply
##' it to \code{fit$model.response}.
##' @title States of a multi-state model
##' @param object Object of class \code{prodlim} or \code{Hist} .
##' @param ... not used
##' @return A character vector with the states of the model.
##' @author Thomas A. Gerds
#' @export
getStates <- function(object,...){
  UseMethod("getStates",object)
}
#' @export 
getStates.Hist <- function(object,...){
  attr(object,"states")
}

#' @export 
getStates.prodlim <- function(object,...){
  attr(object$model.response,"states")
}

