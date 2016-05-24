#' Hessian computation for fitted distance detection function model parameters
#'
#' Computes hessian to be used for variance-covariance matrix.  The hessian is
#' the outer product of the vector of first partials (see pg 62 of Buckland et
#' al 2002).
#'
#' @param ddfobj distance sampling object
#' @param misc.options width-transect width (W); int.range-integration range
#'   for observations; showit-0 to 3 controls level of iteration printing;
#'   integral.numeric-if TRUE integral is computed numerically rather
#'   than analytically
#' @return variance-covariance matrix of parameters in the detection function
#' @note This is an internal function used by \code{\link{ddf.ds}} to fit
#'   distance sampling detection functions.  It is not intended for the user to
#'   invoke this function but it is documented here for completeness.
#' @author Jeff Laake
#' @seealso \code{\link{flnl}},\code{\link{flpt.lnl}},\code{\link{ddf.ds}}
#' @references Buckland et al. 2002
#' @keywords utility
flt.var <- function(ddfobj, misc.options){

  fpar1 <- getpar(ddfobj)
  fpar <- fpar1
  parmat <- NULL

  #   Compute first partial (numerically) of log(f(y)) for each observation
  #   for each parameter and store in parmat (n by length(fpar))
  for (i in 1:length(fpar)){
    deltap <- .0001*fpar1[i]
    if(deltap==0) deltap <- 0.0001
    fpar[i] <- fpar1[i]- deltap
    x1 <- -flpt.lnl(fpar, ddfobj,misc.options)
    fpar[i] <- fpar1[i]+deltap
    x2 <- -flpt.lnl(fpar, ddfobj,misc.options)
    parmat <- cbind(parmat,(x2-x1)/(2*deltap))
  }

  # Compute varmat using first partial approach (pg 62 of Buckland et al 2002)
  varmat <- matrix(0,ncol=length(fpar1),nrow=length(fpar1))

  for(i in 1:length(fpar1)){
    for(j in 1:length(fpar1)){
      varmat[i,j] <- sum(parmat[,i]*parmat[,j])
    }
  }

  return(varmat)
}
