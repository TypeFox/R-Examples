getSIMoNeNet <-
function (DEGeneExpr,NEdges=NA,ClusterMethod="both",AddAnnotations=F,MartDataset="hsapiens_gene_ensembl")
{
	optionValue<-getOption("stringsAsFactors")
	options(stringsAsFactors=F)
	if (class(DEGeneExpr)!="DEGeneExpr")
	{
		stop("Network not of class DEGeneExpr")
	}
	if (!(ClusterMethod %in% c("both","yes","no")))
	{
		stop("Not valid ClusterMethod (both, yes or no)")
	}
	if (!is.numeric(NEdges) & !is.na(NEdges) & !(NEdges %in% c("BIC","AIC")))
	{
		stop("Wrong choice for NEdges")
	}
	
	if (requireNamespace("simone",quietly=TRUE))
	{
	opts=simone::setOptions(verbose=F)
	SIMoNeResults=simone::simone(DEGeneExpr$DataExpression,control=opts) 
	maxBICnedges=SIMoNeResults$n.edges[which.max(SIMoNeResults$BIC)] 
	maxAICnedges=SIMoNeResults$n.edges[which.max(SIMoNeResults$AIC)] 
	TakeEdges=NEdges
	if(is.na(NEdges)){TakeEdges=floor(mean(c(maxBICnedges,maxAICnedges)))} 
	SIMoNeGraph=simone::getNetwork(SIMoNeResults,TakeEdges) 
	
	if (ClusterMethod %in% c("both","yes"))
	{
		opts=simone::setOptions(clusters.crit=TakeEdges,verbose=F) 
		SIMoNeResultsCluster=simone::simone(DEGeneExpr$DataExpression,clustering=T,control=opts)
		maxBICnedges=SIMoNeResultsCluster$n.edges[which.max(SIMoNeResultsCluster$BIC)]
		maxAICnedges=SIMoNeResultsCluster$n.edges[which.max(SIMoNeResultsCluster$AIC)]
		if(is.na(NEdges)){TakeEdges=floor(mean(c(maxBICnedges,maxAICnedges)))}
		SIMoNeGraphCluster=simone::getNetwork(SIMoNeResultsCluster,TakeEdges)
		if (ClusterMethod=="both")
		{
			CancelPos<-which(SIMoNeGraphCluster$A==0)
			SIMoNeGraph$A[CancelPos]<-0
			SIMoNeGraph$A[CancelPos]<-0
			SIMoNeGraph$name<-paste("Common between global and cluster (",length(which(SIMoNeGraph$A==1))/2,")",sep="")
		}
		if (ClusterMethod=="yes") {SIMoNeGraph<-SIMoNeGraphCluster}
	}
	}
	else
	{
		stop("simone package must be installed to use this function")
	}
	
	GenesAnnotations <- NULL
	if (AddAnnotations)
	{
		cat("Adding annotations...\n")
		if (requireNamespace("biomaRt",quietly=TRUE)) {ensembl <- biomaRt::useMart("ensembl", dataset = MartDataset)} else {stop("biomaRt package must be installed to use this function")}
		GenesAnnotations <- getGenesInformations(rownames(SIMoNeGraph$A),ensembl)
	}
	
	SIMoNeNetwork<-SIMoNeNet(SIMoNeGraph,DEGeneExpr,GenesAnnotations)
	options(stringsAsFactors=optionValue)
	return(SIMoNeNetwork)
}
