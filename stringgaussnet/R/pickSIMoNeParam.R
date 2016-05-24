pickSIMoNeParam <-
function (DEGeneExpr,ClusterMethod=F,NEdges=NA)
{
	if (requireNamespace("simone",quietly=TRUE)) {opts=simone::setOptions(verbose=F)} else {stop("simone package must be installed to use this function")}
	if (requireNamespace("simone",quietly=TRUE)) {SIMoNeResults=simone::simone(DEGeneExpr$DataExpression,control=opts)} else {stop("simone package must be installed to use this function")}
	if (ClusterMethod)
	{
		maxBICnedges=SIMoNeResults$n.edges[which.max(SIMoNeResults$BIC)] 
		maxAICnedges=SIMoNeResults$n.edges[which.max(SIMoNeResults$AIC)] 
		TakeEdges=NEdges
		if(is.na(NEdges)){TakeEdges=floor(mean(c(maxBICnedges,maxAICnedges)))}
		if (requireNamespace("simone",quietly=TRUE)) {opts=simone::setOptions(clusters.crit=TakeEdges,verbose=F)} else {stop("simone package must be installed to use this function")}
		if (requireNamespace("simone",quietly=TRUE)) {SIMoNeResults=simone::simone(DEGeneExpr$DataExpression,control=opts)} else {stop("simone package must be installed to use this function")}
		plot(SIMoNeResults)
	}
	if (!ClusterMethod)
	{
		plot(SIMoNeResults)
	}
	
	line <- readline("Press return to close the plot...")
	dev.off()
}
