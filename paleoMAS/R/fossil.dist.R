fossil.dist <-
function(modern,fossil,method="canberra")
{
library(vegan)
{
	modern<-modern$predicted
	names.m<-colnames(modern)
	names.f<-colnames(fossil)
	intersect(names.m,names.f)->names.a
	modern.int<-modern[,names.a]
	fossil.int<-fossil[,names.a]
	distances<-matrix(nrow=nrow(fossil),ncol=nrow(modern))
	for(i in 1:nrow(fossil))
	{
		rbind(fossil.int[i,],modern.int)->base
		as.matrix(vegdist(base,method=method))->dist.base
		dist.base[-1,1]->distances[i,]
	}
	rbind(as.vector(modern[,1]),distances)->distan
}
return(distan)
}

