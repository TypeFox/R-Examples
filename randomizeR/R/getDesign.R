# --------------------------------------------
# Generic function for the Design
# --------------------------------------------

#' Design of a randomization procedure
#' 
#' Generates a \code{character} vector which specifies the used randomization method
#' 
#' @param obj object of the class \code{randSeq} or \code{randPar}.
#'
#' @name getDesign
NULL


#' @rdname getDesign
#' 
#' @export
setGeneric("getDesign", function(obj) standardGeneric("getDesign"))