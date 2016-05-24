#' Calculates projections of the Data
#' 
#' Calculates the projections of the data object onto the directions from an
#' underlying \code{ebblab} object.
#' 
#' The default projection direction is the direction with the best objective
#' criterion.
#' 
#' @name fitted.epplab
#' @aliases fitted.epplab fitted-method fitted,epplab-method
#' @docType methods
#' @param object Object of class \code{epplab}.
#' @param which Onto which direction should the new data be projected.
#' @param ... Additional parameters
#' @return A matrix, having in each column the projection onto the direction of
#' a certain run, and in each row the projected value.
#' @author Daniel Fischer
#' @keywords methods print
#' @examples
#' 
#' library(tourr)
#' data(olive)
#' res <- EPPlab(olive[,3:10],PPalg="PSO",PPindex="KurtosisMin",n.simu=10, maxiter=20)
#' 
#' # Projection to the best direction
#' fitted(res)
#' 
#' # Projection to the 1,2,5 best directions:
#' fitted(res,which=c(1,2,5))
#' 
#' @export
`fitted.epplab` <- function(object,which=1,...){
  which <- which[which<=length(object$PPindexVal)]
  
  PPscores <- object$x %*% object$PPdir[,which,drop=FALSE]

  PPscores
} 
