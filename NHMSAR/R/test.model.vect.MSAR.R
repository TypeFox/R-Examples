test.model.vect.MSAR = function(data,simu,lag=NULL){
	
	N.samples = dim(data)[2]
	N.s = N.samples
	N.sim = dim(simu)[2]
	Bsim = N.sim/N.samples
	
	T = dim(data)[1]
	d = dim(data)[3]
		X.mat = matrix(data,T*N.samples,d)
		X.mat = scale(X.mat,center=TRUE,scale=FALSE)
		Cov.data = t(X.mat)%*%X.mat/dim(X.mat)[1]
		X.mat = matrix(simu,T*N.sim,d)
		X.mat = scale(X.mat,center=TRUE,scale=FALSE)
		Cov.sim = t(X.mat)%*%X.mat/dim(X.mat)[1]
		dd = sum((Cov.data-Cov.sim)^2)
		dd.sim = NULL
		for (k in 1:Bsim) {
			 X.mat = matrix(simu[,((k-1)*N.s+1):(k*N.s),],T*N.samples,d)
			 X.mat = scale(X.mat,center=TRUE,scale=FALSE)
			 Cov.tmp = t(X.mat)%*%X.mat/dim(X.mat)[1]
			 dd.sim[k] = sum((Cov.tmp-Cov.sim)^2)
		}	
	Cvect = list()
	Cvect$dd = dd
	Cvect$q.dd = quantile(dd.sim,probs=c(.05,.95))
	r = rank(c(dd,dd.sim))[1]
	Cvect$p.value = 1-r/(Bsim+1)
	
	return(list(Cvect=Cvect))
}