#------------------------------------------------------------#

### returns the estimated log likelihood of tuning data, for linear learning.

loglik.svm = function(x.train,y.train,x.test,y.test,lambda=0.5,Inum=Inum,
	kernel = kernel, kparam = kparam) 
{
   prob = prob.svm(x.train,y.train,x.test,lambda=lambda,Inum=Inum,
				kernel = kernel, kparam = kparam)
   prob[which(y.test==-1)]=1-prob[which(y.test==-1)]
   prob[prob<=0] = 1e-3
   loglik= - sum( log(prob) ) / nrow(x.test)
   return(loglik)
} ### minus loglikelihood for prob estimation 


#------------------------------------------------------------#