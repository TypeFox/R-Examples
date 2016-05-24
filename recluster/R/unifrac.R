unifrac<-function(data, phylo) {
	util.pd <- beta.pd.utils(data,phylo)
	betadiv <- betajac.pd(util.pd)
	results=NULL 
	results <- dist.mat(data,betadiv[,1])
	names(results)<- colnames(betadiv)[1]
	results
}

