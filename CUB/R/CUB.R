#' @title Main function for CUB models 
#' @description Main function to estimate and validate a CUB model for explaining uncertainty 
#' and feeling for given ratings, with or without covariates and shelter effect.
#' @aliases CUB
#' @usage CUB(ordinal,m=get('m',envir=.GlobalEnv),
#' Y=0,W=0,shelter=0,maxiter,toler,makeplot=TRUE)
#' @param ordinal Vector of ordinal responses
#' @param m Number of ordinal categories (if omitted, it will be assigned to the number of categories
#'  specified in the global environment)
#' @param Y Matrix of covariates for explaining the uncertainty component. If omitted (default), no covariate 
#' is included in the model
#' @param W Matrix of covariates for explaining the feeling component. If omitted (default), no covariate
#' is included in the model
#' @param shelter Category corresponding to the shelter choice. If omitted (default), no shelter effect
#' is included in the model
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm 
#' (default: maxiter=500) 
#' @param toler Fixed error tolerance for final estimates (default: toler = 1e-4)
#' @param makeplot Logical: if TRUE (default), the algorithm returns a graphical plot comparing fitted 
#' probabilities and observed relative frequencies for CUB models without covariates. If only one 
#' explicative dichotomous variable is included in the model for the feeling or the uncertainty components, 
#' then the function returns a graphical plot comparing the distributions of the responses conditioned to
#'  the value of the covariate
#' @export CUB
#' @return An object of the class "CUB" is a list containing the following results: 
#' \item{estimates}{Maximum likelihood estimates: \eqn{(\pi, \xi)}}
#' \item{loglik}{Log-likelihood function at the final estimates}
#' \item{varmat}{Variance-covariance matrix of final estimates}
#' \item{niter}{Number of executed iterations}
#' \item{BIC}{BIC index for the estimated model}
#' @details This is the main function for CUB models, which calls for the corresponding functions whenever 
#' covariates or shelter effect are specified. It performs maximum likelihood estimation via the E-M algorithm 
#' for CUB models and extensions. The optimization procedure is run via "optim".\cr
#' It is possible to fit data with CUB models, with or without covariates 
#' for the parameters of the mixture model, and CUB models with shelter effect with no covariate included 
#' in the model. The program also checks if the estimated variance-covariance matrix is positive definite: 
#' if not, it prints a warning message and returns a matrix and related results with NA entries.
#' @references 
#' Piccolo D. and D'Elia A. (2008). A new approach for modelling consumers' preferences, \emph{Food Quality and Preference},
#' \bold{18}, 247--259 \cr
#' Iannario M. and Piccolo D. (2012). CUB models: Statistical methods and empirical evidence, in: 
#' Kenett R. S. and Salini S. (eds.), \emph{Modern Analysis of Customer Surveys: with applications using R}, 
#' J. Wiley and Sons, Chichester, 231--258\cr
#' Iannario M. (2012). Modelling \emph{shelter} choices in a class of mixture models for ordinal responses,  
#' \emph{Statistical Methods and Applications}, \bold{21}, 1--22 \cr
#' Iannario M. and Piccolo D. (2014). Inference for CUB models: a program in R, \emph{Statistica & Applicazioni}, 
#' \bold{XII} n.2, 177--204
#' @seealso \code{\link{probcub00}}, \code{\link{probcubp0}}, \code{\link{probcub0q}}, \code{\link{probcubpq}},
#' \code{\link{probcubshe1}}, \code{\link{loglikCUB}}, \code{\link{varmatCUB}} 
#' @keywords models
#' @examples 
#' \donttest{
#' data(relgoods)
#' m<-10
#' ordinal<-relgoods[,40] 
#' model<-CUB(ordinal)     # Equivalent calls: CUB(ordinal, m) or CUB(ordinal,m=10) 
#'      # if m has not been previously declared
#' estpar<-model$estimates  # Estimated parameter vector (pai,csi)
#' maxlik<-model$loglik     # Log-likelihood function at ML estimates
#' vmat<-model$varmat
#' nniter<-model$niter
#' BICCUB<-model$BIC
#' ################
#' ## CUB model
#' data(univer)
#' m<-7
#' officeho<-univer[,10]
#' model<-CUB(officeho)
#' BICcub<-model$BIC
#' ################
#' ## CUB model with covariate for uncertainty
#' data(relgoods)
#' m<-10
#' ordinal<-relgoods[,26] 
#' gender<-relgoods[,7]
#' data<-na.omit(cbind(ordinal,gender))
#' modelcovpai<-CUB(data[,1],Y=data[,2])
#' BICcov<-modelcovpai$BIC
################
#' ## CUB model with covariate for feeling
#' data(univer)
#' m<-7
#' ordinal<-univer[,12]
#' freqserv<-univer[,2]
#' modelcovcsi<-CUB(ordinal,W=freqserv)
#' ##################
#' ## CUB model with covariates for both components
#' data(univer)
#' m<-7
#' gender<-univer[,4]
#' lage<-log(univer[,3])-mean(log(univer[,3]))
#' ordinal<-univer[,12]
#' maxiter<-500; toler<-1e-6;
#' model<-CUB(ordinal,Y=gender,W=lage) # Makeplot is ignored
#' param<-model$estimates
#' bet<-param[1:2]      # ML estimates of coefficients for uncertainty covariate
#' gama<-param[3:4]     # ML estimates of coefficients for feeling covariate
#' }


CUB <-function(ordinal,m=get('m',envir=.GlobalEnv),Y=0,W=0,shelter=0,maxiter,toler,makeplot=TRUE){
  
  if (missing(maxiter)){
    maxiter <- 500
  }
  if (missing(toler)){
    toler <- 1e-4
  }
  
  
  ry<-NROW(Y);   rw<-NROW(W); shelter<-as.numeric(shelter)
  if(shelter!=0){
    if (ry==1 & rw==1){
      cubshe(m,ordinal,shelter,maxiter,toler,makeplot)
    } else {
      cat("CUB model with shelter effect available only with no covariates")
    }
  } else{
    if(ry==1 & rw==1) cub00(m,ordinal,maxiter,toler,makeplot)
    else{
      if(ry!=1 & rw==1) cubp0(m,ordinal,Y,maxiter,toler,makeplot)
      else{
        if(ry==1 & rw!=1) cub0q(m,ordinal,W,maxiter,toler,makeplot)
        else{
          if(ry!=1 & rw!=1) cubpq(m,ordinal,Y,W,maxiter,toler,makeplot)
          else cat("Wrong variables specification")
        }
      }                            
    }
  }
}
