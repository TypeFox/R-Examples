#' @title Variance-covariance matrix for CUBE models
#' @aliases varmatCUBE
#' @description Compute the variance-covariance matrix of parameter estimates for CUBE models when no covariate
#' is specified, or when covariates are included for all the three parameters.
#' @usage varmatCUBE(ordinal,m,param,Y=0,W=0,Z=0,expinform=FALSE)
#' @export varmatCUBE
#' @param ordinal Vector of ordinal responses
#' @param m Number of ordinal categories
#' @param param Vector of parameters for the specified CUBE model
#' @param Y Matrix of selected covariates to explain the uncertainty component (default: no covariate is included 
#' in the model)
#' @param W Matrix of selected covariates to explain the feeling component (default: no covariate is included 
#' in the model)
#' @param Z Matrix of selected covariates to explain the overdispersion component (default: no covariate is included 
#' in the model)
#' @param expinform Logical: if TRUE  and no covariate is included in the model, the function returns
#'  the expected variance-covariance matrix (default is FALSE: the function returns the observed 
#'  variance-covariance matrix)
#' @details The function checks if the variance-covariance matrix is positive-definite: if not, 
#' it returns a warning message and produces a matrix with NA entries.
#' @seealso  \code{\link{CUBE}}, \code{\link{loglikCUBE}}
#' @keywords htest
#' @references Iannario, M. (2014). Modelling Uncertainty and Overdispersion in Ordinal Data, 
#' \emph{Communications in Statistics - Theory and Methods}, \bold{43}, 771--786 \cr
#' Piccolo, D. (2014). Inferential issues on CUBE models with covariates, 
#' \emph{Communications in Statistics - Theory and Methods}, \bold{44}, DOI: 10.1080/03610926.2013.821487
#' @examples
#' m<-7; n<-500;
#' pai<-0.83; csi<-0.19; phi<-0.045;
#' ordinal<-simcube(n,m,pai,csi,phi)
#' param<-c(pai,csi,phi)
#' varmat<-varmatCUBE(ordinal,m,param)
#' ##########################
#' ### Including covariates
#' \donttest{
#' data(relgoods)
#' m<-10
#' ordinal<-relgoods[,37]
#' age<-2014-relgoods[,4]
#' lage<-log(age)-mean(log(age))
#' nona<-na.omit(cbind(ordinal,lage))
#' ordinal<-nona[,1]
#' Y<-W<-Z<-nona[,2]
#' estbet<-c(0.18, 1.03); estgama<-c(-0.6, -0.3); estalpha<-c(-2.3,0.92);
#' param<-c(estbet,estgama,estalpha);
#' varmat<-varmatCUBE(ordinal, m, param, Y=Y, W=W, Z=Z, expinform=TRUE)
#' }


varmatCUBE <-
function(ordinal,m,param,Y=0,W=0,Z=0,expinform=FALSE){
   ry<-NROW(Y); rw<-NROW(W); rz<-NROW(Z);


  if(ry==1 & rw==1 & rz==1) {
    pai<-param[1]; csi<-param[2]; phi<-param[3];
    if(expinform==FALSE){
     freq<-tabulate(ordinal,nbins=m)
     varmat<-varcovcubeobs(m,pai,csi,phi,freq)
     } else{
       n<-length(ordinal)
        varmat<-varcovcubeexp(m,pai,csi,phi,n)
      }
  } else {
    if(ry>1 & rz>1 & rw >1){
      ncy<-NCOL(Y)
     ncw<-NCOL(W)
     estbet<-param[1:(ncy+1)]; estgama<-param[(ncy+2):(ncy+ncw+2)]; estalpha<-param[(ncy+ncw+3):length(param)];
     varmat<-varcovcubecov(m,ordinal,Y,W,Z,estbet,estgama,estalpha) 
     } else {
       cat("CUBE models not available for this variables specification")
     }
  }
return(varmat)

}
