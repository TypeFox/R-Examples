computeSimilarities <-
function (DEGeneExpr)
{
	Correlations=PValues=matrix(NA,nrow=ncol(DEGeneExpr$DataExpression),ncol=ncol(DEGeneExpr$DataExpression)) 
	rownames(Correlations)=colnames(Correlations)=rownames(PValues)=colnames(PValues)=colnames(DEGeneExpr$DataExpression) 
	
	for (i in 1:nrow(Correlations))
	{
		for (j in 1:nrow(Correlations))
		{
			gene1=rownames(Correlations)[i]
			gene2=colnames(Correlations)[j]
			TestSpearman=pspearman::spearman.test(DEGeneExpr$DataExpression[,gene2],DEGeneExpr$DataExpression[,gene1],approximation="AS89") 
			Correlations[gene1,gene2]=TestSpearman$estimate
			PValues[gene1,gene2]=TestSpearman$p.value
		}
	}
	diag(Correlations)<-0 
	diag(PValues)<-1
	
	Similarities<-(1+Correlations)/2 
	
	Values<-list(Correlations=Correlations,PValues=PValues,Similarities=Similarities)
	return(Values)
}
