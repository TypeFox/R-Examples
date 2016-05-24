matrixBC.BIC <-
function(x,k,r,lambda,alpha=0.2,beta=0.2,nstart=20,Sigma.init=NULL,Delta.init=NULL){
	x<-x-mean(x)
	BIC<-rep(NA,length(lambda))
    nonzero<-rep(NA,length(lambda))

	for(i in 1:length(lambda)){
		mclustering<-matrixBC(x,k,r,lambda=lambda[i],nstart=nstart,alpha=alpha,beta=beta,Sigma.init=Sigma.init,Delta.init=Delta.init)
		BIC[i]<-CalculateBICMatrix(x,mclustering)
		nonzero[i]<-sum(mclustering$Mus!=0)

	}
	return(list(lambda=lambda[which(BIC==min(BIC))[1]],BIC=BIC,nonzeromus=nonzero))
}
