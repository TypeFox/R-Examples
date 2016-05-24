sparseBC.BIC <-
function(x,k,r,lambda){
	x<-x-mean(x)
	BIC<-rep(NA,length(lambda))
	nonzero<-rep(NA,length(lambda))
	for(i in 1:length(lambda)){
		bires<-sparseBC(x,k,r,lambda=lambda[i])
		BIC[i]<-CalculateBIC(x,bires)
		nonzero[i]<-sum(bires$Mus!=0)
	}
	return(list(lambda=lambda[which(BIC==min(BIC))[1]],BIC=BIC,nonzeromus=nonzero))
}
