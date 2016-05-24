# gradients of the likelihood for line transects
flt.gr<-function(fpar,mix.terms,width,x,ftype="hn",z=NULL,zdim=0,showit=0,pt=FALSE,EM=FALSE){
   ### First juggle the parameters.
   swpar<-switchpars(fpar,mix.terms,z=z,zdim=zdim)
   fpar<-swpar$fpar
   z<-swpar$z
   zdim<-swpar$zdim
   got.pars<-getpars(fpar,mix.terms,zdim,z)
   key.scale<-got.pars$key.scale
   key.shape<-got.pars$key.shape
   mix.prop<-c(got.pars$mix.prop)

   # use the correct multipliers and integration functions when using
   # point transects
   if(pt){
      mult<-2*pi*x$distance
      intfcn<-integrate.hn.pt
   }else{
      mult<-1
      intfcn<-integrate.hn
   }

   # shorthand for number of samples...
   n<-length(x$distance)
   # detection function eval at every distance
   g<-detfct(x$distance,fpar,mix.terms,zdim,z)
   # overall mu/nu
   int<-mu.calc(fpar,mix.terms,width,z,zdim,pt)
   if(length(int)==1){
      int<-rep(int,n)
   }

   ##################
   # beta derivatives

   intercepts<-c()

   # intercept terms
   for(j in 1:mix.terms){

      if(mix.terms>1){
         mix.propj<-mix.prop[j]
      }else{
         mix.propj<-mix.prop
      }

      if(is.list(z)|all(zdim>0)){
         if(length(z)==1){
            zz<-z[[1]]
         }else{
            zz<-z[[j]]
         }
         keysc<-key.scale[j,]
      }else{
         keysc<-key.scale[j]
         zz<-matrix(1,length(x$distance),1)
      }

      intj<-intfcn(keysc,width)
      if(length(intj)==1){
         intj<-rep(intj,n)
      }

      # first fb, detection function derivative
      # same for point and line transects
      fb<-(x$distance^2)/(keysc^2)*
            keyfct.hn(x$distance,keysc)

      # second bit -- different for lt and pt
      if(pt){
         sb<-4*pi*keysc^2*(1-keyfct.hn(width,keysc))-
               2*pi*width^2*keyfct.hn(width,keysc)
      }else{
         sb<-(intj-width*keyfct.hn(width,keysc))
      }

      intercepts<-c(intercepts,mix.propj*sum(fb/g-sb/int))
   }
   dbeta<-intercepts

   # common terms
   # number of common pars is # pars - # intercepts - # alphas
   K<-length(fpar)-mix.terms-(mix.terms-1)
   if(K>0){
      # loop over covariates
      for(k in 1:K){

         fb<-sb<-rep(0,n)

         for(j in 1:mix.terms){
            if(is.list(z)|all(zdim>0)){
               if(length(z)==1){
                  zz<-z[[1]]
               }else{
                  zz<-z[[j]]
               }
               keysc<-key.scale[j,]
            }else{
               keysc<-key.scale[j]
               zz<-matrix(1,length(x$distance),1)
            }

            intj<-intfcn(keysc,width)

            fb<-fb+mix.prop[j]*zz[,k+1]*
                 (x$distance/keysc)^2*
                 keyfct.hn(x$distance,keysc)

            # second bit -- different for lt and pt
            if(pt){
               sb<-sb+mix.prop[j]*2*(zz[,k+1]*intj-
                        zz[,k+1]*pi*width^2*keyfct.hn(width,keysc))
            }else{
               sb<-sb+mix.prop[j]*zz[,k+1]*
                      (intj-width*keyfct.hn(width,keysc))
            }
         }
         dbeta<-c(dbeta,sum(fb/g-sb/int))
      }
   }


   ##############################
   # derivatives wrt mixture prop
   if(mix.terms>1){
      dpi<-rep(0,(mix.terms-1))
      alphas<-fpar[(length(fpar)-mix.terms+2):length(fpar)]


      # a function for the (A_j - A_{j-1}) terms
      Af<-function(j,js,alphas){
         esum1<-esum2<-0
         if((j>= js)){
            esum1<-dgamma(sum(exp(alphas[1:j])),3,2)
         }
         if((j-1)>=js){
            esum2<-dgamma(sum(exp(alphas[1:(j-1)])),3,2)
         }
         return(exp(alphas[js])*(esum1-esum2))
      }


      if(is.list(z)|all(zdim>0)){
         keyscJ<-key.scale[mix.terms,]
      }else{
         keyscJ<-key.scale[mix.terms]
      }

      gJ<-keyfct.hn(x$distance,keyscJ)
      intJ<-intfcn(keyscJ,width)

      # js is the alpha we're finding the derivative wrt
      # THINK js = j*
      for(js in 1:(mix.terms-1)){
         # first bit -- the detection function bit
         fb<-rep(0,n)
         # second bit -- the mu bit
         sb<-rep(0,n)
         Jterm<-0

         for(j in 1:(mix.terms-1)){
            if(is.list(z)|all(zdim>0)){
               keysc<-key.scale[j,]
            }else{
               keysc<-key.scale[j]
            }

            # what is A_j - A_{j-1}
            Afj<-Af(j,js,alphas)

            fb<-fb+Afj*(keyfct.hn(x$distance,keysc)-gJ)
            mudiff<-(intfcn(keysc,width)-intJ)

            if(length(mudiff)==1){
               mudiff<-rep(mudiff,n)
            }
            sb<-sb+Afj*mudiff
         }

         dpi[js]<-sum(fb/g-sb/int)
      }
   }else{
      dpi<-c()
   }
   
   #################
   # return vector
   dpar<-c(dbeta,dpi)

   if(showit>2)
      cat("Derivative pars=",dpar,"\n")

   return(dpar)
}
