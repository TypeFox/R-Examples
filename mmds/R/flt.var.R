"flt.var" <- function(fpar, x,width, mix.terms,z=NULL,zdim=0,pt=FALSE)
#   flt.var - computes hessian for v-c matrix (see pg 62 of Buckland et al 2002)
#   Value:  variance-covariance matrix of parameters
#   Functions Used:  flt
{
#
 fpar1<-fpar
  parmat=NULL

eps<-0.0001

#   Compute first partial (numerically) of log(f(y)) for each observation 
#   for each parameter and store in parmat (n by length(fpar))
  for (i in 1:length(fpar)){
    fpar[i] <- fpar1[i]-eps*fpar1[i]
	x1=log(rowSums(eval.pdf(fpar,x,width,mix.terms,0,"hn",z,zdim,pt)))
    fpar[i] <- fpar1[i]+eps*fpar1[i]
	x2=log(rowSums(eval.pdf(fpar,x,width,mix.terms,0,"hn",z,zdim,pt)))
    parmat=cbind(parmat,(x2-x1)/(2*eps*fpar1[i]))
  }
#
#    Compute varmat using first partial approach (pg 62 of Buckland et al 2002)
#
  varmat=matrix(0,ncol=length(fpar1),nrow=length(fpar1))

  for(i in 1:length(fpar1))
    for(j in 1:length(fpar1))
      varmat[i,j]=sum(parmat[,i]*parmat[,j]) 
  
  return(varmat)
}
