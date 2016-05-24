integrate.pdf<-function(x,width,pars,mix.terms,z=NULL,zdim=0,pt=FALSE){
   # integrate the detection function -- ie. calculate overall mu

   # grab some pars
   got.pars<-getpars(pars,mix.terms,zdim,z)
   key.scale<-got.pars$key.scale
   key.shape<-got.pars$key.shape
   mix.prop<-got.pars$mix.prop

   if(pt){
      intfcn<-integrate.hn.pt
   }else{
      intfcn<-integrate.hn
   }

   #work out prob det for each mixture
   #if(mix.terms>1){
      if(is.list(z)|all(zdim>0)){

         p<-matrix(0,mix.terms,length(x$distance))
         for (i in 1:mix.terms) {

            if(is.list(z)){
               keysc<-key.scale[i,]
            }else{
               keysc<-key.scale[i]
            }

           p[i,]<-intfcn(keysc,width)
         }
      }else{
         p<-numeric(mix.terms)
         for (i in 1:mix.terms) {

            if(is.list(z)){
               keysc<-key.scale[i,]
            }else{
               keysc<-key.scale[i]
            }

           p[i]<-intfcn(keysc,width)
         }
      }
      

      #work out proportion of each mixture class in the population
      p.pop<-mix.prop/p
      p.pop<-p.pop/sum(p.pop)
   #}else{
   #   p.pop<-mix.prop
   #}

   res<-numeric(length(x$distance))
   # covariates...
   if(is.list(z)|all(zdim>0)){
      # storage
   
      for(j in 1:mix.terms){
         res<-res+(p.pop[j,]/intfcn(key.scale[j,],width))*intfcn(key.scale[j,],x$distance)
      }
   # no covariates
   }else{
      # storage
      for(j in 1:mix.terms){
         res<-res+(p.pop[j]/intfcn(key.scale[j],width))*intfcn(key.scale[j],x$distance)
      }
   }
   return(res)
}
