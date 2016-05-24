unifract<-function(data, phylo) {
	util.pd <- beta.pd.utils(data,phylo)
	betadiv <- betajac.pd(util.pd)
	results=NULL 
	results <- dist.mat(data,betadiv[,2])
	names(results)<- colnames(betadiv)[2]
	results
}

