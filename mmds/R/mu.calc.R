mu.calc<-function(pars,mix.terms,width,z=NULL,zdim=0,pt=FALSE){
   # calculate mu across mixtures 
   # returns 1 number for non-covariate mixtures
   # n-vector for covariate mixtures

   gp<-getpars(pars,mix.terms,zdim,z)
   sigma<-gp$key.scale
   pis<-gp$mix.prop

   if(pt){
      intfcn<-integrate.hn.pt
   }else{
      intfcn<-integrate.hn
   }

   if(mix.terms>1){
      if(is.null(z)){
         mu<-sum(pis*intfcn(sigma,width))
      }else{
         mu<-apply(sigma,1,intfcn,width)
         mu<-mu%*%matrix(pis,length(pis),1)
      }
   }else{
      if(is.null(z)){
         mu<-intfcn(sigma,width)
      }else{
         mu<-apply(sigma,1,intfcn,width)
      }
   }
   return(mu)
}
