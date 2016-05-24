# --------------------------------------------
# Generic function for the lineup
# --------------------------------------------

#' Lineup of a united formation
#' 
#' Generates a \code{numeric} vector which specifies the used united lineup
#' 
#' @param obj object of the class \code{formation}.
#' 
#' @return \code{vector} of the used lineup
#'
#' @name getLineup
NULL

#' @rdname getLineup
#' 
#' @export
setGeneric("getLineup", function(obj) standardGeneric("getLineup"))