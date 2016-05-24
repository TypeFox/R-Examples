#' Calculates projections for a new Data Object
#' 
#' Calculates the projections of a new data object onto the directions from an
#' existing \code{ebblab} object.
#' 
#' The default projection direction is the direction with the best objective
#' criterion. In case that no data is given to the function, the fitted scores
#' for the original data will be returned.
#' 
#' @name predict.epplab
#' @aliases predict.epplab predict-method predict,epplab-method
#' @docType methods
#' @param object Object of class \code{epplab}.
#' @param which Onto which direction should the new data be projected.
#' @param data The new data object
#' @param ... Additional parameters
#' @return A matrix having in each column the projection onto the direction of
#' a certain run and in each row the projected value.
#' @author Daniel Fischer
#' @keywords methods print
#' @examples
#' 
#' library(tourr)
#' data(olive)
#' res <- EPPlab(olive[,3:10], PPalg="PSO", PPindex="KurtosisMin", n.simu=10, maxiter=20)
#' 
#' newData <- matrix(rnorm(80), ncol=8)
#' 
#' # Projection on the best direction
#' predict(res, data=newData)
#' 
#' # Projection on the best 3 directions
#' predict(res, which=1:3, data=newData)
#' 
#' # Similar with function fitted() when no data is given:
#' predict(res)
#' fitted(res)
#' 
#' @export
`predict.epplab` <- function(object,which=1, data=NULL,...){
  
  which <- which[which<=length(object$PPindexVal)]
  if(is.null(data)==T) data <- object$x

  # Center the data with the mean vector of the object:
  data <- as.matrix(sweep(data,2,object$center,"-"))

  PPscores <- data %*% object$PPdir[,which]

  PPscores
} 
