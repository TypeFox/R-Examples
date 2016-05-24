#########################################################################################################################################################
##########################################################  Generic functions   #########################################################################
#########################################################################################################################################################

##############################################################
#' Continuing estimation.
#'
#' Generic function, that given the results of an estimation procedure (ex. MCMC or maximum likelihood optimization)
#' continues the procedure for some more iterations.
#'
#' @param object An object representing the results of an estimation procedure which we wish to continue.
#'               For example it might represents an MCMC chain.
#'               
#' @param ... Additional arguments that might be needed to continue the estimation procedure.
#'
#' @return An object of the same class as \code{object}, where the results of the estimation have been updated.
#' 
#' @seealso For examples, see \code{\link{smcmc-class}}.
#' 
#' @export
#' @docType methods
#' @rdname continue-generic

setGeneric('continue', function(object, ...) standardGeneric("continue"))


