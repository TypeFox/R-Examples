namepars<-function(fpar,mix.terms,zdim=0,z=NULL){
   ### name the parameter values
   # Args:
   #  fpar          parameters from optimisation
   #  mix.terms     number of mixture components
   #  zdim          dimension of the z matrix/matrices
   #  z             covariate matrix (matrices)


   # now split for the various cases

   ### In the non-covariate case
   if(is.null(zdim)){
      beta.labels<-paste("beta_",1:mix.terms,sep="")
   }else if(all(zdim==0)){
      beta.labels<-paste("beta_",1:mix.terms,sep="")
   ### multiple covariate case with 1-point mixture
   }else if(mix.terms==1){
      
      z<-z[[1]]
      # when we are calculating P_a z is not a matrix, so be careful
      if(is.null(dim(z))){
         ncol<-length(z)
      }else{
         ncol<-dim(z)[2]
      }
      beta.labels<-paste("beta_",0:(ncol-1),sep="")
   ### multiple covariate mixtures
   }else{
      
      ### multiple covariate mixtures - all pars vary
      if(length(zdim)>1){
         i<-1
         # loop over the zs
         for (zmat in z){
            beta.labels<-paste("beta_",i,0:(zdim[i]-1),sep="")
            i<-i+1
         }
      ### common pars, vary intercept
      }else{
         # should be supplied c(intercept pars, other pars)
         # want to pick different intercept par each time and 
         # the same common pars each time

         beta.labels<-paste("beta_",1:mix.terms,0,sep="")
   
         beta.labels<-c(beta.labels,paste("beta_",1:(zdim-1),sep=""))
      }
   }

   ## mixture proportion names alpha_1... alpha_{J-1}
   # don't do anything if we only have a 1 point mixture
   if(mix.terms != 1){
      mix.prop.labels<-paste("alpha_",1:(mix.terms-1),sep="")
      beta.labels<-c(beta.labels,mix.prop.labels)
   }

   #return the labels
   return(beta.labels)
}
