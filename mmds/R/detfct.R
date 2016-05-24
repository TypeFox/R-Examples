detfct<-function(x,pars,mix.terms,zdim=0,z=NULL){
   # evaluate the total detection function

   gp<-getpars(pars,mix.terms,zdim,z)
   sigmas<-gp$key.scale
   pis<-gp$mix.prop

   res<-rep(0,length(x))

   for(j in 1:mix.terms){

      if(is.null(z)){
         keysc<-sigmas[j]
      }else{
         keysc<-sigmas[j,]
      }

      res<-res+pis[j]*keyfct.hn(x,keysc)
   }
   res
}
