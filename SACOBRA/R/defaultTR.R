#defaultTR.r
#'
#' Default settings for the trust-region part of COBRA.
#' 
#' Sets default values for the trust-region part \code{cobra$TRlist} of COBRA.  
#' With the call \code{\link{setOpts}(myTR,defaultTR())} it is possible to extend a partial list 
#' \code{myTR} to a list containing all \code{TR}-elements (the missing ones are taken from 
#' \code{defaultTR()}).

#'
#'  @return a list with the following elements 
#'    \describe{
#'      \item{radiMin}{Minimum fraction of the width of the search space to be used as radius of the trust region [c(0:1)]}
#'      \item{radiMax}{Maximum fraction of the width of the search space to be used as radius of the trust region [c(0:1)]}
#'      \item{radiInit}{Initial radius of trust region}
#'}
#'


defaultTR<-function(){
  tr<-list(radiMin=0.01,
           radiMax=0.5,
           radiInit=0.05)
  return(tr)  
}

