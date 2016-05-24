#' Root Expected Square Difference
#' @param x a column vector of scores on which the rsd is conditioned
#' @param o a column vector of equated scores based on the overall population (aligned with elements in x)
#' @param g a column vector of equated scores based on a single subpopulation (aligned with elements in x)
#' @param f a column vector of relative frequency associated with each raw score (can be based on either overall population or a subpopulation) (aligned with elements in x)
#' @param s a scalar representing the standard deviation of x for any (sub)population of interest (e.g., synthetic population) (default is 1, which leads to calculation of the unstandardized resd)
#' @author Anne Corinne Huggins-Manley 
#' @references 
#' \itemize{
#' \item{Yang, W.L. (2004). Sensitivity of linkings between AP multiple-choice scores and composite scores to geographical region: An illustration of checking for population invariance. Journal of Educational Measurement, 41, 33-41.}
#' }
#' @seealso \code{\link{rsd}}
#' @examples
#' #Unstandardized RESD for subpopulation 1 in the example data set, ex.data
#' resd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,3],f=ex.data[,8])
#' 
#' #Unstandardized RESD for subpopulation 5 in the example data set, ex.data
#' resd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,7],f=ex.data[,8])
#' 
#' #Standardized RESD for subpopulation 5 in the example data set, ex.data
#' resd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,7],f=ex.data[,8],s=4.2)
#' 
#' @return root expected square difference
#' @export

resd<- function(x,o,g,f,s){
  
  if(missing(s))
    s <- 1
  
  resdX<-(o-g)^2
  resdY<-f*resdX
  resdZ<-(sqrt(sum(resdY)))/s  
  return(resdZ)
}