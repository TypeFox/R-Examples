llnhlogit=function(theta,choice,lnprices,Xexpend) 
{
#   function to evaluate non-homothetic logit likelihood
#       choice is a n x 1 vector with indicator of choice (1,...,m)
#	lnprices is n x m array of log-prices faced
#       Xexpend is n x d array of variables predicting expenditure
#
#       non-homothetic model specifies ln(psi_i(u))= alpha_i - exp(k_i)u
#
#       structure of theta vector:
#           alpha  (m x 1)
#           k      (m x 1)
#           gamma  (k x 1)   expenditure function coefficients
#           tau   scaling of v
#	    
   m=ncol(lnprices)
   n=length(choice)
   d=ncol(Xexpend)
   alpha=theta[1:m]
   k=theta[(m+1):(2*m)]
   gamma=theta[(2*m+1):(2*m+d)]
   tau=theta[length(theta)]
   iotam=c(rep(1,m))
   c1=as.vector(Xexpend%*%gamma)%x%iotam-as.vector(t(lnprices))+alpha
   c2=c(rep(exp(k),n))
   u=callroot(c1,c2,.0000001,20)
   v=alpha - u*exp(k)-as.vector(t(lnprices))
   vmat=matrix(v,ncol=m,byrow=TRUE)
   vmat=tau*vmat
   ind=seq(1,n)
   vchosen=vmat[cbind(ind,choice)]
   lnprob=vchosen-log((exp(vmat))%*%iotam)
   return(sum(lnprob))
}
