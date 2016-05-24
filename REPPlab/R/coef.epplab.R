#' Extracts the Directions of an Epplab Object
#' 
#' Extracts the found directions of an \code{epplab} object.
#' 
#' The coef function extracts the directions found from the EPPlab call.
#' 
#' @name coef.epplab
#' @aliases coef.epplab coef-method coef,epplab-method
#' @docType methods
#' @param object Object of class \code{epplab}.
#' @param which Specifies which directions are extracted.
#' @param ... Additional parameters.
#' @author Daniel Fischer
#' @keywords methods print
#' @examples
#' 
#' library(tourr)
#' data(olive)
#' res <- EPPlab(olive[,3:10],PPalg="PSO",PPindex="KurtosisMin",n.simu=10, maxiter=20)
#' coef(res)
#' 
#' @export
`coef.epplab` <- function(object,which=1:ncol(object$PPdir),...){
 object$PPdir[,which,drop=FALSE]
} 
