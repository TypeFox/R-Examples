#' @title Main function for CUBE models
#' @description Main function to estimate and validate a CUBE model for given ratings, 
#' explaining uncertainty, feeling and overdispersion.
#' @aliases CUBE
#' @usage CUBE(ordinal, m=get('m',envir=.GlobalEnv), Y = 0, W = 0, Z = 0, starting, maxiter,
#' toler, makeplot = TRUE, expinform=FALSE)
#' @param ordinal Vector of ordinal responses
#' @param m Number of ordinal categories (if omitted, it will be assigned to the number of categories
#'  specified in the global environment)
#' @param Y Matrix of covariates for explaining the uncertainty component. If omitted (default), no covariate 
#' is included in the model
#' @param W Matrix of covariates for explaining the feeling component. If omitted (default), no covariate
#' is included in the model
#' @param Z Matrix of covariates for explaining the overdispersion component. If omitted (default), no covariate
#' is included in the model
#' @param starting  Initial parameters estimates to start the optimization algorithm. If missing,
#'  the function calls specific routines computing the best initial estimates available
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
#' (default: maxiter=1000) 
#' @param toler Fixed error tolerance for final estimates (default: toler = 1e-6,
#' except for the model including covariates for all the three parameters, in which case toler=1e-2)
#' @param makeplot Logical: if TRUE (default) and no covariate is included in the model, the function returns
#'  a graphical plot comparing fitted probabilities and observed relative frequencies
#' @param expinform Logical: if TRUE  and no covariate is included in the model, the function returns 
#'  the expected variance-covariance matrix (default is FALSE)
#' @export CUBE
#' @return An object of the class "CUBE" is a list containing the following results: 
#' \item{estimates}{Maximum likelihood estimates: \eqn{(\pi, \xi, \phi)}}
#' \item{loglik}{Log-likelihood function at the final estimates}
#' \item{varmat}{Variance-covariance matrix of final estimates}
#' \item{niter}{Number of executed iterations}
#' \item{BIC}{BIC index for the estimated model}
#' @details It is the main function for CUBE models, calling for the corresponding functions whenever
#'  covariates are specified: it is possible to select covariates for explaining all the three parameters
#'   or only the feeling component. \cr
#' The program also checks if the estimated variance-covariance matrix is positive definite: if not,
#'  it prints a warning message and returns a matrix and related results with NA entries.
#'  The optimization procedure is run via "optim". If covariates are included only for feeling,
#' the variance-covariance matrix is computed as the inverse of the returned numerically differentiated
#'  Hessian matrix (option: hessian=TRUE as argument for "optim"), and the estimation procedure is not
#'  iterative, so a NULL result for $niter is produced.
#'  If the estimated variance-covariance matrix is not positive definite, the function returns a 
#'  warning message and produces a matrix with NA entries.
#' @references 
#' Iannario, M. (2014). Modelling Uncertainty and Overdispersion in Ordinal Data, 
#' \emph{Communications in Statistics - Theory and Methods}, \bold{43}, 771--786 \cr
#' Piccolo, D. (2014). Inferential issues on CUBE models with covariates,
#'  \emph{Communications in Statistics - Theory and Methods}, \bold{44}, DOI: 10.1080/03610926.2013.821487 \cr
#' Iannario, M. (2015). Detecting latent components in ordinal data with overdispersion by means
#'  of a mixture distribution, \emph{Quality & Quantity}, \bold{49}, 977--987
#' @seealso \code{\link{probcube}}, \code{\link{loglikCUBE}}, \code{\link{loglikcuben}},  \code{\link{inibestcube}},
#'  \code{\link{inibestcubecsi}}, \code{\link{inibestcubecov}},
#' \code{\link{varmatCUBE}}
#' @keywords models
#' @examples 
#' \donttest{
#' data(relgoods)
#' m<-10
#' ordinal<-na.omit(relgoods[,37])  
#' model<-CUBE(ordinal,starting=c(0.1,0.1,0.1))  
#' model$estimates        # Final ML estimates
#' model$loglik           # Maximum value of the log-likelihood function
#' model$varmat         
#' model$niter
#' model$BIC
#' ######################## 
#' ordinal<-relgoods[,40]
#' cov<-relgoods[,2]
#' nona<-na.omit(cbind(ordinal,cov))
#' modelcovcsi<-CUBE(nona[,1],W=nona[,2])
#' modelcov<-CUBE(nona[,1],Y=nona[,2],W=nona[,2], Z=nona[,2])
#' modelcov$BIC
#' modelcovcsi$BIC
#' #######################################
#' data(univer)
#' m<-7
#' ordinal<-univer[,8]
#' starting<-inibestcube(m,ordinal)
#' model<-CUBE(ordinal,starting=starting)
#' }

CUBE <-function(ordinal,m=get('m',envir=.GlobalEnv),Y=0,W=0,Z=0,starting,maxiter,toler,
               makeplot=TRUE,expinform=FALSE){#default FAlse for expinform
  
  ry<-NROW(Y); rw<-NROW(W); rz<-NROW(Z);
  
  
  if(missing(starting)){
    if(ry==1 & rw==1 & rz==1) {
      starting<-inibestcube(m,ordinal)
    }else if(ry==1 & rz==1 & rw >1){
      initial<-inibestcube(m,ordinal)
      starting<-inibestcubecsi(m,ordinal,W,initial,maxiter=500,toler=1e-6)
    } else {
      starting<-inibestcubecov(m,ordinal,Y,W,Z)
    }
  }
  
  
  if(missing(maxiter)){
    maxiter<-1000
  }
  if(missing(toler)){
    toler<-1e-6
  }
  if(missing(expinform)){
    expinform<-FALSE
  }
  
  if(ry==1 & rw==1 & rz==1) {
    cube000(m,ordinal,starting,maxiter,toler,makeplot,expinform)
  } else if (ry==1 & rz==1 & rw >1){
    cubecsi(m,ordinal,W,starting,maxiter,toler)  
  } else if(ry>1 & rz>1 & rw >1){
    cubecov(m,ordinal,Y,W,Z,starting,maxiter,toler=1e-2)
  } else {
    cat("CUBE models not available for this variables specification")
  }
}
