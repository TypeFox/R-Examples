"flt"<-function(fpar, x,width, mix.terms,showit=0, ftype="hn",z=NULL,zdim=0,pt=FALSE,EM=FALSE){
   # flt - a function to optimize to find the MLE 
   # Arguments
   #  fpar        parameter values
   #  x           distances
   #  width       truncation width
   #  mix.terms   number of mixture terms
   #  showit      level of debugging 0=none up to 3
   #  ftype       function type (only "hn" at the moment)
   #  z           covariate matrix
   #  zdim        how many covariates
   # Returns
   #  lnl         log-likelihood   

   # evaluate the pdf at each datum
   res<-eval.pdf(fpar,x,width,mix.terms,showit=showit, 
                    ftype=ftype,z=z,zdim=zdim,pt=pt,EM=EM)

   # eval.pdf returns a (n obs)x(mix.terms) matrix
   res<-rowSums(res)

   # Calculate the sum of the log likelihood
   lnl <- sum(log(res),na.rm=TRUE)

   # debug output -- log likelihood and pars
   if(showit>2){
     cat("lt lnl = ", lnl,   "\n")    
     cat("par = ", fpar,"\n")
   }
  
   return(lnl)
}
