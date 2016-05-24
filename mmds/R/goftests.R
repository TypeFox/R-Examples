"goftests"<-function(x,width,pars,mix.terms,z=NULL,zdim=0,pt=FALSE){
	# Kolmogorv-Smirnov  and Cramer-von Mises tests for 
   # mixture detection functions.
	# Mostly hacked together from mrds implementation.
	# Arguments
	#  x					vector of distances
	#  width				truncation width
	#  pars				detection function parameters
	#  mix.terms		the number of mixture terms
	# Return -- a list
   # $ks$
	#     p			Kolmogorov-Smirnov test p-value
	#     D_n		Kolmogorov-Smirnov test statistic
   # $cvm$
	#     p			Cramer-von Mises test p-value
	#     W  		Cramer-von Mises test statistic

   if(pt){
      intfcn<-integrate.hn.pt
   }else{
      intfcn<-integrate.hn
   }
   got.pars<-getpars(pars,mix.terms,zdim,z)
   key.scale<-got.pars$key.scale
   key.shape<-got.pars$key.shape
   mix.prop<-got.pars$mix.prop

   cdf.eval<-rep(0,length(x$distance))

   for(i in 1:length(x$distance)){

      for(j in 1:mix.terms){
         if(is.list(z)|all(zdim>0)){
            keysc<-key.scale[j,i]
            keysh<-key.shape[j,i]
         }else{
            keysc<-key.scale[j]
            keysh<-key.shape[j]
         }

         cdf.eval[i]<-cdf.eval[i]+mix.prop[j]*intfcn(keysc,x$distance[i])
      }
   }

   cdf.eval<-cdf.eval/mu.calc(pars,mix.terms,width,z,zdim,pt)

   # sort by function evaluation
   cdfvalues<-sort(cdf.eval)
   n<-length(cdfvalues)

   # EDF indicator function
   edf.ind <- function(x,a,lt){
      if(lt)
         length(a[a<x])
      else 
         length(a[a<=x])
   }

   # evaluate the EDF 
   lower.edf<-(unlist(sapply(cdfvalues,edf.ind,a=cdfvalues,lt=TRUE)))/n
   upper.edf<-(unlist(sapply(cdfvalues,edf.ind,a=cdfvalues,lt=FALSE)))/n
   ################  

   ########## KS ##############
   # calculate the test statistic
   Dn<-max(max(abs(lower.edf-cdfvalues)),max(abs(upper.edf-cdfvalues)))

   # if the max returned -Inf (something bad happened)
   # then just return NA
   if(is.infinite(Dn)){
      Dn<-NA
      ks.p<-NA
   }else{
	   # This is taken directly from pks in mrds
      diff<-1
      ks.p<-exp(-2*n*Dn^2)
      i<-1
      while(abs(diff)>.0000001){
         i<-i+1
         diff<-(-1)^(i-1)*exp(-2*n*(i*Dn)^2)
         ks.p<-ks.p+diff
      }
   }

   ########## C-vM ############
   W<-1/(12*n) + sum((cdfvalues - ((1:n)-.5)/n)^2)

   log.eps <- log(1e-05)
   y <- matrix(0, nrow = 4, ncol = length(W))
   for (k in 0:3) {
      z <- gamma(k+0.5)*sqrt(4*k+1)/(gamma(k+1)*pi^(3/2)*sqrt(W))
      u <- (4*k+1)^2/(16*W)
      y[k+1, ] <- ifelse(u >-log.eps, 0, z * exp(-u)*besselK(x=u,nu=1/4))
   }
   cvm.p<-1-apply(y, 2, sum)

   # bundle some objects with nice classes so the print methods below work...
   ks<-list(p=2*ks.p,Dn=Dn)
   cvm<-list(p=cvm.p,W=W)

   class(ks)<-"ds.mixture.ks"
   class(cvm)<-"ds.mixture.cvm"

   # return the test statistics and p-values
   return(list(ks=ks,cvm=cvm))
}


# plotting methods
print.ds.mixture.ks<-function(x,...){
   cat("\nKolmogorov-Smirnov test\n")
   cat("Dn =",x$Dn," p-value =",x$p,"\n")
}
print.ds.mixture.cvm<-function(x,...){
   cat("\nCramer-von Mises test\n")
   cat("W =",x$W," p-value =",x$p,"\n")
}
