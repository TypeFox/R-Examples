#' Reparametrize mixture proportions
#'
#' Re-parameterise the mixture proportions so that when there is more than a 
#' 2-point mixture, the proportions sum to 1.
#'
#' @param thetas Mixture proportions in their parametrisation for optimization.
#' @return a vector of parameters 
#'
#' @export
#'
#' @section Notes:
#' See Miller and Thomas for information on exactly how these are calculated. 
#' Thanks go to David Borchers for proposing the trick.
#' 
#' @references
#' Miller, D.L. and L. Thomas (in prep.). Mixture model distance sampling detection functions.
#'
#' @seealso inv.reparam.pi
#'
#' @author David L. Miller
#' @examples
#' library(mmds)
#' reparam.pi(inv.reparam.pi(0.3))
#' reparam.pi(inv.reparam.pi(c(0.3,0.4,0.1),lastpar=TRUE))
reparam.pi<-function(thetas){
   # re-parameterise the pis 

   # function to use for the transform
   F<-function(x) pgamma(x,3,2) # gamma cdf
   #F<-function(x){exp(x)/(1+exp(x))} # logistic

   # so we don't have to constrain to be >0
   ethetas<-exp(thetas)

   mix.props<-rep(0,length(thetas))

   for(i in 1:length(ethetas)){

      lastbit<-F(sum(ethetas[1:(i-1)]))
      if(i==1) lastbit<-0

      mix.props[i]<-F(sum(ethetas[1:i]))-lastbit
   }

   # return the pis
   return(c(mix.props,1-sum(mix.props)))
}
