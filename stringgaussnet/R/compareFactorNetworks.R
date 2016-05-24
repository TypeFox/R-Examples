compareFactorNetworks <-
function (Networks,Colors=rainbow(length(Networks)),interactiveMode=T)
{
	if (class(Networks)!="FactorNetworks"){stop("Bad class used instead of FactorNetworks.")}
	MergedEdges<-mergeFactorEdges(Networks)
	
	interactivePrompt<-"Press return for next plot..."
	
	Levels<-unique(MergedEdges$Level)
	
	Connecs<-vector()
	for (Level in Levels)
	{
		LevelEdges<-MergedEdges[which(MergedEdges$Level==Level),]
		Connecs[Level]<-mean(table(c(LevelEdges$node1,LevelEdges$node2)))
	}
	
	Rhos<-abs(MergedEdges$Rho)
	PValues<--log(MergedEdges$P.Value,10)
	maxpval<-99
	finitepos<-which(is.finite(PValues))
	if(length(finitepos)>0){maxpval<-max(PValues[finitepos])}
	infinitepos<-which(!is.finite(PValues))
	if(length(infinitepos)>0){PValues[infinitepos]<-maxpval}
	
	boxplot(Rhos~MergedEdges$Level,col=Colors,main="Absolute values of rho in different networks",ylab="abs(Rho)")
	
	if(interactiveMode){line <- readline(interactivePrompt)}
	boxplot(PValues~MergedEdges$Level,col=Colors,main="Spearman p-values in different networks",ylab="-log10(P-values)")
	
	if(interactiveMode){line <- readline(interactivePrompt)}
	barplot(table(MergedEdges$Level),col=Colors,main="Number of edges in different networks",ylab="number of edges")
	
	if(interactiveMode){line <- readline(interactivePrompt)}
	barplot(Connecs,col=Colors,main="Connectivities in different networks",ylab="mean(connectivity)")
	
	if(interactiveMode){line <- readline("Press return to close the plot...");dev.off()}
}
