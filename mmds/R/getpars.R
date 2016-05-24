#' Grab parameter values
#'
#' Extract parameter values and create a named list.
#'
#' @param fpar The \code{$par} element \code{\link{ds.mixture}} object.
#' @param mix.terms Number of mixture components.
#' @param zdim Number of covariates.
#' @param z Covariate matrix.
#' @return a named list with elements \code{$key.scale} (giving the key scales)
#'     and \code{$mix.prop} giving the mixture proportions.
#'
#' @author David L. Miller
getpars<-function(fpar,mix.terms,zdim=0,z=NULL){
   # grab the mixture proportions
   if(mix.terms != 1){
      mix.prop.tmp<-fpar[(length(fpar)-mix.terms+2):length(fpar)]

      # do the clever re-parameterisation here
      mix.prop<-reparam.pi(mix.prop.tmp)

      #mix.prop<-c(mix.prop.tmp,1-sum(mix.prop.tmp))
      fpar<-fpar[-((length(fpar)-mix.terms+2):length(fpar))]
   }else{
      mix.prop<-as.vector(1)
   }

   # now split for the various cases

   ### In the non-covariate case
   if(is.null(zdim)){
      # get all the key scales
      key.scale<-exp(fpar)
   }else if(all(zdim==0)){
      # get all the key scales
      key.scale<-exp(fpar)

   ### multiple covariate case with 1-point mixture
   }else if(mix.terms==1){
      
      z<-z[[1]]

      # when we are calculating P_a z is not a matrix, so be careful
      if(is.null(dim(z))){
         ncol<-length(z)
      }else{
         ncol<-dim(z)[1]
      }

      # calculating the key scale is a bit of a pain.
      # end up returning a matrix with # cols = mix terms
      #                                # rows = # data
      key.scale<-matrix(NA,nrow=1,ncol=ncol)
      key.scale[1,]<-scalevalue(fpar,z)
   ### multiple covariate mixtures
   }else{
      # structure to hold the parameters
      key.scale<-matrix(NA,nrow=mix.terms,ncol=dim(z[[1]])[1])
      
      ### multiple covariate mixtures - all pars vary
      if(length(zdim)>1){
         i<-1
         # loop over the zs
         for (zmat in z){
            # select the parameters from fpar
            these.pars<-fpar[1:zdim[i]]
            key.scale[i,]<-scalevalue(these.pars,zmat)

            # get rid of what we've already used
            fpar<-fpar[(zdim[i]+1):length(fpar)]

            i<-i+1
         }
      ### common pars, vary intercept
      }else{
         # should be supplied c(intercept pars, other pars)
         # want to pick different intercept par each time and 
         # the same common pars each time

         # grab the intercepts, put the other terms in opars
         intercepts<-fpar[1:mix.terms]
         opars<-fpar[(mix.terms+1):length(fpar)]

         for(j in 1:mix.terms){
            key.scale[j,]<-scalevalue(c(intercepts[j],opars),z[[1]])
         }
      }
   }

   # Return the result
   return(list(key.scale=key.scale,mix.prop=mix.prop))
}
