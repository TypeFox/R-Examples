#' @title Variance-covariance matrix for CUB models
#' @aliases varmatCUB
#' @description Compute the variance-covariance matrix of parameter estimates for CUB models with or without
#' covariates for the feeling and the overdispersion parameter, and for extended CUB models with shelter effect.
#' @usage varmatCUB(ordinal,m,param,Y=0,W=0,shelter=0)
#' @export varmatCUB
#' @param ordinal Vector of ordinal responses
#' @param m Number of ordinal categories
#' @param param Vector of parameters for the specified CUB model
#' @param Y Matrix of selected covariates to explain the uncertainty component (default: no covariate is included 
#' in the model)
#' @param W Matrix of selected covariates to explain the feeling component (default: no covariate is included 
#' in the model)
#' @param shelter Category corresponding to the shelter choice (default: no shelter effect is included in the 
#' model)
#' @details The function checks if the variance-covariance matrix is positive-definite: if not, 
#' it returns a warning message and produces a matrix with NA entries.
#' @seealso  \code{\link{CUB}}, \code{\link{loglikCUB}}
#' @keywords htest
#' @references Piccolo D. (2006). Observed Information Matrix for MUB Models, 
#' \emph{Quaderni di Statistica}, \bold{8}, 33--78 \cr
#' Iannario, M. (2012). Modelling shelter choices in ordinal data surveys. 
#' \emph{Statistical Modelling and Applications}, \bold{21}, 1--22
#' @examples
#' data(univer)
#' m<-7
#' ### CUB model with no covariate
#' ordinal<-univer[,12]
#' pai<-0.87; csi<-0.17; 
#' param<-c(pai,csi)
#' varmat<-varmatCUB(ordinal,m,param)
#' #######################
#' ### and With covariates for feeling
#' data(univer)
#' m<-7
#' ordinal<-univer[,9]
#' pai<-0.86; gama<-c(-1.94, -0.17);
#' param<-c(pai,gama)
#' W<-univer[,4]           
#' varmat<-varmatCUB(ordinal,m,param, W=W)
#' #######################
#' ### CUB model with uncertainty covariates
#' data(relgoods)
#' m<-10
#' ordinal<-relgoods[,29]
#' gender<-relgoods[,2]
#' data<-na.omit(cbind(ordinal,gender))
#' ordinal<-data[,1]
#' Y<-data[,2]
#' bet<-c(-0.811,0.93); csi<-0.202;
#' varmat<-varmatCUB(ordinal,m,param=c(bet,csi),Y=Y)
#' #######################
#' ### and with covariates for both parameters
#' data(relgoods)
#' m<-10
#' gender<-relgoods[,2]
#' smoking<-relgoods[,12]
#' ordinal<-relgoods[,40]
#' nona<-na.omit(cbind(ordinal,gender,smoking))
#' ordinal<-nona[,1]
#' gender<-nona[,2]
#' smoking<-nona[,3]
#' gama<-c(-0.55, -0.43); bet<-c(-0.45, -0.48);
#' varmat<-varmatCUB(ordinal,m,param=c(bet,gama),Y=gender,W=smoking)
#' #######################
#' ### Variance-covariance for a CUB model with shelter
#' m<-8; n<-300;
#' pai1<-0.5; pai2<-0.3; csi<-0.4;
#' shelter<-6;
#' pr<-probcubshe1(m,pai1,pai2,csi,shelter)
#' ordinal<-sample(1:m,n,prob=pr,replace=TRUE)
#' param<-c(pai1,pai2,csi)
#' varmat<-varmatCUB(ordinal,m,param,shelter=shelter)

varmatCUB <-
function(ordinal,m,param,Y=0,W=0,shelter=0){

    ry<-NROW(Y);   rw<-NROW(W); shelter<-as.numeric(shelter)
  
  if(shelter!=0){
    if (ry==1 & rw==1){
      pai1<-param[1]
      pai2<-param[2]
      csi<-param[3]
      n<-length(ordinal)
      varmat<-varcovcubshe(m,pai1,pai2,csi,shelter,n)
    } else {
      cat("CUB model with shelter effect available only with no covariates")
    }
  }else{
    if(ry==1 & rw==1) {
      pai<-param[1]
      csi<-param[2]
      varmat<-varcovcub00(m,ordinal,pai,csi)
    }    else{
      if(ry!=1 & rw==1) {
       ncy<-NCOL(Y)
       bet<-param[1:(ncy+1)]
       csi<-param[length(param)]
       varmat<-varcovcubp0(m,ordinal,Y,bet,csi)
       } else {
        if(ry==1 & rw!=1) {
        pai<-param[1]
        gama<-param[2:length(param)]
        varmat<-varcovcub0q(m,ordinal,W,pai,gama)
       } else{
          if(ry!=1 & rw!=1) {
            ncy<-NCOL(Y)
            bet<-param[1:(ncy+1)]
          gama<-param[(ncy+2):length(param)]
          varmat<-varcovcubpq(m,ordinal,Y,W,bet,gama)
        } else {
           cat("Wrong variables specification")
        }
      }                            
    }
  }
  }
  return(varmat)
}
