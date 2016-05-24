#' Plots a Scatterplot Matrix for an epplab Object
#' 
#' Plots a scatterplot matrix of fitted scores of an \code{epplab} object.
#' 
#' The option \code{which} can restrict the output to certain simulation runs.
#' In case of many simulations, this might improve the readability.
#' 
#' @name pairs.epplab
#' @aliases pairs.epplab pairs-method pairs,epplab-method
#' @docType methods
#' @param x Object of class \code{epplab}.
#' @param which Which simulation runs should be taken into account.
#' @param ... Graphical parameters, see also par().
#' @author Daniel Fischer
#' @keywords methods hplot
#' @examples
#' 
#' library(tourr)
#' data(olive)
#' res <- EPPlab(olive[,3:10],PPalg="PSO",PPindex="KurtosisMin",n.simu=10, maxiter=20)
#' pairs(res)
#' 
#' @export
`pairs.epplab` <- function(x,which=1:10,...){
  
  which <- which[which<=length(x$PPindexVal)]
  x.fitted <- fitted(x,which=which)
  
  pairs(x.fitted,...)

  invisible()
} 
