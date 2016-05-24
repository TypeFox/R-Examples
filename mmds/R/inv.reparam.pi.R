#' Reparametrize mixture proportions
#'
#' Re-parameterise the mixture proportions so that when there is more than a 
#' 2-point mixture, the proportions sum to 1.
#'
#' @param mix.prop A set of mixture proportions.
#' @param lastpar Is the last parameter provided, i.e. does sum(mix.prop)=1?
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
#' @seealso reparam.pi
#'
#' @author David L. Miller
#' @examples
#' library(mmds)
#' reparam.pi(inv.reparam.pi(0.3))
#' reparam.pi(inv.reparam.pi(c(0.3,0.4,0.1),lastpar=TRUE))
inv.reparam.pi<-function(mix.prop,lastpar=FALSE){
   # re-parameterise the pis so that when there is more
   # than a 2 point mixture
   # INVERSE VERSION!
   # pis->thetas

   # functions to use for the transforms
   F<-function(x) pgamma(x,3,2)
   Finv<-function(x) qgamma(x,3,2)
   #F<-function(x){exp(x)/(1+exp(x))}
   #Finv<-function(x){log(x/(1-x))}

   thetas<-c()

   if(lastpar){
      looper<-length(mix.prop)
   }else{
      looper<-(length(mix.prop)-1)
   }

   for(i in 1:looper){
      if(i==1){
         theta<-log(Finv(mix.prop[i]))
      }else{
         theta<-log(
                  Finv(mix.prop[i]+F(sum(exp(thetas))))-
                  sum(exp(thetas)))
      }
      thetas<-c(thetas,theta)
   }

   # return the thetas
   return(thetas)
}
