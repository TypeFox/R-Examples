# --------------------------------------------
# Generic function for simulating a red card
# --------------------------------------------

#' Simulate red card(s)
#' 
#' Simulates red card(s) in the united and returns the adjusted lineup.
#' 
#' @param obj object of the class \code{formation}
#' @param lineup lineup of the corresponding object \code{obj}
#' 
#' @return \code{vector} of the adjusted lineup for the red card(s)
#'
#' @name simRedCard
NULL

#' @rdname simRedCard
#' 
#' @export
setGeneric("simRedCard", function(obj, lineup) standardGeneric("simRedCard"))