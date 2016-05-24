#------------------------------------------------------------#

### returns the prob estimation for linear learning.

prob.svm = function(x.train,y.train,x.test,lambda ,Inum ,
	kernel , kparam  ) 
{
	 
	I = seq(0,1,length=Inum+2)[2:(Inum+1)]
	
	obj = pi.path(lambda = lambda*nrow(x.train), x.train, y.train, kernel  = kernel, kparam = kparam, eps = 1e-5, ridge=1e-10)
	
	tempprob = predict_pipathresults(obj,x.test,pi=I)
	res<-I[1]/2+apply(tempprob$predicted.y+1,1,sum)/2/(Inum+1)
	return(res)
	
} ### svm learning, prob calculation for test data