getWGCNANet <-
function (DEGeneExpr,SoftThreshold=8,AThreshold=0.85,AddAnnotations=F,MartDataset="hsapiens_gene_ensembl")
{
	optionValue<-getOption("stringsAsFactors")
	options(stringsAsFactors=F)
	CorrelationVals<-computeSimilarities(DEGeneExpr)
	Correlations<-CorrelationVals$Correlations
	PValues<-CorrelationVals$PValues
	Similarities<-CorrelationVals$Similarities
	rm(CorrelationVals)
	AdjacencyScores=1/(1+exp(-SoftThreshold*(Similarities-0.5))) 
	AdjacencyScores[which(AdjacencyScores<AThreshold & AdjacencyScores>(1-AThreshold))]<-0
	AdjacencyScores[lower.tri(AdjacencyScores,diag=T)]<-0
	
	GenesAnnotations <- NULL
	if (AddAnnotations)
	{
		cat("Adding annotations...\n")
		if (requireNamespace("biomaRt",quietly=TRUE)) {ensembl <- biomaRt::useMart("ensembl", dataset = MartDataset)} else {stop("biomaRt package must be installed to use this function")}
		GenesAnnotations <- getGenesInformations(rownames(AdjacencyScores),ensembl)
	}
	
	WGCNANetwork<-WGCNANet(AdjacencyScores,SoftThreshold,AThreshold,Correlations,PValues,DEGeneExpr,GenesAnnotations)
	options(stringsAsFactors=optionValue)
	return(WGCNANetwork)
}
