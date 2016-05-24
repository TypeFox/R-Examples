#' @title Log-likelihood function for CUB models
#' @aliases loglikCUB
#' @description  Compute the log-likelihood function for CUB models with or without covariates to 
#' explain the feeling and uncertainty components, or for extended CUB models with shelter effect.
#' @usage loglikCUB(ordinal,m,param,Y=0,W=0,shelter=0)
#' @export loglikCUB
#' @param ordinal Vector of ordinal responses
#' @param m Number of ordinal categories
#' @param param Vector of parameters for the specified CUB model
#' @param Y Matrix of selected covariates to explain the uncertainty component (default: no covariate is included 
#' in the model)
#' @param W Matrix of selected covariates to explain the feeling component (default: no covariate is included 
#' in the model)
#' @param shelter Category corresponding to the shelter choice (default: no shelter effect is included in the 
#' model)
#' @details If no covariate is included in the model, then "param" has the form \eqn{(\pi,\xi)}. More generally, 
#' it has the form \eqn{(\beta,\gamma)} where, respectively, \eqn{\beta} and \eqn{\gamma} are the vectors of 
#' coefficients explaining the uncertainty and the feeling components, with length NCOL(Y)+1 and
#'  NCOL(W)+1 to account for an intercept term in the first entry.
#' @seealso  \code{\link{CUB}}, \code{\link{cubforsim}}
#' @keywords htest
#' @examples
#' ## Log-likelihood of a CUB model with no covariate
#' m<-9; n<-300;
#' pai<-0.6; csi<-0.4;
#' ordinal<-simcub(n,m,pai,csi)
#' param<-c(pai,csi)
#' loglikcub<-loglikCUB(ordinal,m,param)
#' ##################################
#' ## Log-likelihood of a CUB model with covariate for uncertainty
#' \donttest{
#' data(relgoods)
#' m<-10; 
#' ordinal<-relgoods[,29]
#' gender<-relgoods[,2]
#' data<-na.omit(cbind(ordinal,gender))
#' ordinal<-data[,1]; Y<-data[,2];
#' bbet<-c(-0.81,0.93); ccsi<-0.2;
#' param<-c(bbet,ccsi)
#' loglikcubp0<-loglikCUB(ordinal,m, param, Y=Y)
#' #######################
#' ## Log-likelihood of a CUB model with covariate for feeling
#' data(relgoods)
#' m<-10; 
#' ordinal<-relgoods[,29]
#' gender<-relgoods[,2]
#' data<-na.omit(cbind(ordinal,gender))
#' ordinal<-data[,1]; W<-data[,2];
#' pai<-0.44; gama<- c(-0.91, -0.7);
#' param<-c(pai,gama)
#' loglikcub0q<-loglikCUB(ordinal,m,param,W=W)
#' }
#' #######################
#' ## Log-likelihood of a CUB model with covariates for both parameters
#' data(relgoods)
#' m<-10;
#' smoking<-relgoods[,12]
#' ordinal<-relgoods[,40]
#' gender<-relgoods[,2]
#' nona<-na.omit(cbind(ordinal,gender,smoking))
#' ordinal<-nona[,1]
#' gender<-nona[,2]; smoking<-nona[,3];
#' bet=c(-0.45, -0.48); gama=c(-0.55, -0.43);
#' param<-c(bet,gama)
#' loglikcubpq<-loglikCUB(ordinal,m,param, Y=smoking, W=gender)
#' #################################
#' ### Log-likelihood of a CUB model with shelter effect
#' m<-7; n<-400;
#' pai1<-0.56; pai2<-0.34; csi<-0.16;
#' shelter<-5
#' pr<-probcubshe1(m,pai1,pai2,csi,shelter)
#' ordinal<-sample(1:m,n,prob=pr,replace=TRUE)
#' param<-c(pai1, pai2, csi)
#' loglik<-loglikCUB(ordinal,m,param,shelter=shelter)

loglikCUB <-
function(ordinal,m,param,Y=0,W=0,shelter=0){

  freq<-tabulate(ordinal,nbins=m)
  
  ry<-NROW(Y);   rw<-NROW(W); shelter<-as.numeric(shelter)
  
  if(shelter!=0){
    if (ry==1 & rw==1){
      pai1<-param[1]
      pai2<-param[2]
      csi<-param[3]
      loglik<-loglikcubshe(m,freq,pai1,pai2,csi,shelter)
    } else {
      cat("CUB model with shelter effect available only with no covariates")
    }
  }else{
    if(ry==1 & rw==1) {
      pai<-param[1]
      csi<-param[2]
      loglik<-loglikcub00(m,freq,pai,csi)
    }
    else{
      if(ry!=1 & rw==1) {
        ncy<-NCOL(Y)
       bbet<-param[1:(ncy+1)]
       ccsi<-param[length(param)]
       loglik<-loglikcubp0(m,ordinal,Y,bbet,ccsi) 
       } else {
        if(ry==1 & rw!=1) {
        pai<-param[1]
        gama<-param[2:length(param)]
        loglik<-loglikcub0q(m,ordinal,W,pai,gama)}
        else{
          if(ry!=1 & rw!=1) {
            ncy<-NCOL(Y)
            bet<-param[1:(ncy+1)]
          gama<-param[(ncy+2):length(param)]
          loglik<-loglikcubpq(m,ordinal,Y,W,bet,gama)}
          else cat("Wrong variables specification")
        }
      }                            
    }
  }
  return(loglik)
}
