eval.pdf<-function(fpar,x,width,mix.terms,showit=0,ftype="hn",z=NULL,zdim=0,pt=FALSE,...){
   #### evaluate the pdf

   # get some parameters
   got.pars<-getpars(fpar,mix.terms,zdim,z)
   key.scale<-got.pars$key.scale
   key.shape<-got.pars$key.shape
   mix.prop<-got.pars$mix.prop

   # if the data is from point transects set the key function
   # multiplier to be the distances
   if(pt){
      mult<-2*pi*x$distance
      intfcn<-integrate.hn.pt
   }else{
      mult<-1
      intfcn<-integrate.hn
   }

   # Initialise the result vector
   res<-c()
   # integral -- mu or nu
   int1<-0

   # loop over the mixture parts, summing as we go
   for (j in 1:mix.terms){
      if(is.list(z)|all(zdim>0)){
         keysc<-key.scale[j,]
         keysh<-key.shape[j,]
      }else{
         keysc<-key.scale[j]
         keysh<-key.shape[j]
      }

      # calculate detection function value at the point
      p1<-keyfct(x$distance, keysc,keysh,ftype)

      # calculate the component for this mixture
      # and add a row to the matrix...
      pdf.eval<-mix.prop[j]*p1
      
      int1<-int1+mix.prop[j]*intfcn(keysc,width)

      res<-cbind(res,pdf.eval)
   }

   # done!
   return(mult*res/int1)
}
