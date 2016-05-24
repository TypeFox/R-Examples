# put the pars in the right order
switchpars<-function(fpar,mix.terms,zdim=0,z=NULL){

   # grab the mixture proportions (actually thetas)
   if(mix.terms != 1){
      # all but the last parameter
      mix.prop.tmp<-fpar[(length(fpar)-mix.terms+2):length(fpar)]
      # get the last parameter too
      tmpsum<-sum(reparam.pi(mix.prop.tmp)[1:(mix.terms-1)])
      mix.prop.last<-inv.reparam.pi(c(1-tmpsum,tmpsum))

      mix.prop<-c(mix.prop.tmp,mix.prop.last)

      fpar<-fpar[-((length(fpar)-mix.terms+2):length(fpar))]
   }else{
      mix.prop<-c(1)
      return(list(fpar=fpar,z=z,zdim=zdim)) # if we only have one mixture, exit
   }

   # get the order right
   mix.order<-order(mix.prop)
   mix.prop<-mix.prop[mix.order]

   ########## now split for the various cases

   ### In the non-covariate case
   if(is.null(zdim) || all(zdim==0)){
      # get all the key scales
      fpar<-fpar[mix.order]

   ### multiple covariate mixtures
   }else{

      ### multiple covariate mixtures - all pars vary
      if(length(zdim)>1){
         
         # make a list of the pars for each mixture
         these.pars<-list()
         # loop over the zs
         i<-1
         for (zmat in z){
            # put them into the list in the right order
            these.pars[[mix.order[i]]]<-fpar[1:zdim[i]]
            # get rid of what we've already used
            fpar<-fpar[(zdim[i]+1):length(fpar)]
            i<-i+1
         }
         # make a vector out of it
         fpar<-unlist(these.pars)

         # sort out z and zdim
         zdim<-zdim[mix.order]
         oldz<-z
         j<-1
         for(i in mix.order){
            z[[j]]<-oldz[[i]]
            j<-j+1
         }

      ### common pars, vary intercept
      }else{
         # should be supplied c(intercept pars, other pars)
         # this is easier, just need to change the order
         # of the intercepts
   
         intercepts<-fpar[1:mix.terms]
         the.rest<-fpar[(mix.terms+1):length(fpar)]

         fpar<-c(intercepts[mix.order],the.rest)
      # sort out z and zdim
      if(length(zdim)>1){
         zdim<-zdim[mix.order]
         oldz<-z
         j<-1
         for(i in mix.order){
            z[[j]]<-oldz[[i]]
            j<-j+1
         }
      }
      }
     
   }

   # Return the result
   return(list(fpar=c(fpar,mix.prop[1:(length(mix.prop)-1)]),z=z,zdim=zdim,ord=mix.order))
}
