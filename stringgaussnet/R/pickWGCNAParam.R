pickWGCNAParam <-
function (DEGeneExpr,Alphas=c(c(1:10), seq(from = 12, to=20, by=2)),AThreshold=0.85,interactiveMode=T)
{
	CorrelationVals<-computeSimilarities(DEGeneExpr)
	Correlations<-CorrelationVals$Correlations
	PValues<-CorrelationVals$PValues
	Similarities<-CorrelationVals$Similarities
	rm(CorrelationVals)
	
	Adjacencies<-list()
	for (Alpha in Alphas)
	{
		Adjacencies[[as.character(Alpha)]] = 1/(1+exp(-Alpha*(Similarities-0.5))) 
	}
	
	PossibleColors=rainbow(length(Adjacencies)) 
	names(PossibleColors)=names(Adjacencies)
	
	NEdges=PValuesMax=ConnecMeans=MinRhos=vector()
	for (Alpha in names(Adjacencies))
	{
		Abscissa="Similarity"
		Ordinate="Adjacency"
		Titre="Adjacency computation parameters"
		if (Alpha!=names(Adjacencies)[1]) 
		{
			par(new=T)
			Abscissa=""
			Ordinate=""
			Titre=""
		}
		Adjacency=Adjacencies[[Alpha]]
		plot(c(Adjacency)~c(Similarities),bg=PossibleColors[Alpha],ylim=c(0,1),xlim=c(0,1),pch=21,col="black",main=Titre,xlab=Abscissa,ylab=Ordinate) 
		EdgesPos=which(Adjacency>AThreshold | Adjacency<(1-AThreshold)) 
		NEdges=c(NEdges,length(EdgesPos)/2) 
		EdgesMatrix=matrix(0,ncol=ncol(Adjacency),nrow=nrow(Adjacency)) 
		RhosMatrix=matrix(Inf,ncol=ncol(Correlations),nrow=nrow(Correlations)) 
		colnames(EdgesMatrix)=colnames(RhosMatrix)=colnames(Adjacency) 
		rownames(EdgesMatrix)=rownames(RhosMatrix)=rownames(Adjacency)
		if (length(EdgesPos) > 0) 
		{
			ThresPos=which.min(abs(Adjacency-AThreshold))
			PValuesMax=c(PValuesMax,PValues[ThresPos])
		}
		if (length(EdgesPos) == 0) 
		{
			PValuesMax=c(PValuesMax,1)
		}
		EdgesMatrix[EdgesPos]=Adjacency[EdgesPos] 
		RhosMatrix[EdgesPos]=abs(Correlations[EdgesPos]) 
		Connectivities=apply(EdgesMatrix,2,function(x){return(length(which(x!=0)))}) 
		ConnecMeans=c(ConnecMeans,mean(Connectivities)) 
		MinRhos=c(MinRhos,min(RhosMatrix[which(RhosMatrix!=0)])) 
	}
	names(NEdges)=names(PValuesMax)=names(ConnecMeans)=names(MinRhos)=names(Adjacencies) 
	abline(h=AThreshold,lty=2,lwd=2) 
	abline(h=1-AThreshold,lty=2,lwd=2)
	legend("topleft",legend=names(Adjacencies),pch=21,col="black",pt.bg=PossibleColors,ncol=3,bg="white",title="Alpha values") 
	
	interactivePrompt<-"Press return for next plot..."
	
	if(interactiveMode){line <- readline(interactivePrompt)}
	plot(NEdges~names(NEdges),type="n",xlab="Alpha",ylab="Edges number",main="Edges number as a function of Alpha")
	text(names(NEdges),NEdges,names(NEdges),col="red",cex=2)
	Percents=c(0.05,0.1,0.25) 
	abline(h=Percents*length(Correlations)/2,lty=2,lwd=2)
	
	if(interactiveMode){line <- readline(interactivePrompt)}
	plot(-log(PValuesMax,10)~names(PValuesMax),type="n",xlab="Alpha",ylab="-log10(max(P-Values))",main="Maximal P-Values as a function of Alpha")
	text(names(PValuesMax),-log(PValuesMax,10),names(PValuesMax),col="red",cex=2)
	PValsBreaks=c(0.05,0.01,0.001,0.0001) 
	abline(h=-log(PValsBreaks,10),lty=2,lwd=2)
	
	if(interactiveMode){line <- readline(interactivePrompt)}
	plot(ConnecMeans~names(ConnecMeans),type="n",xlab="Alpha",ylab="mean(Connectivities)",main="Connectivities as a function of Alpha") 
	text(names(ConnecMeans),ConnecMeans,names(ConnecMeans),col="red",cex=2)
	
	if(interactiveMode){line <- readline(interactivePrompt)}
	plot(MinRhos~names(MinRhos),type="n",xlab="Alpha",ylab="min(Rhos)",main="Minimal rhos as a function of Alpha")
	text(names(MinRhos),MinRhos,names(MinRhos),col="red",cex=2)
	RhosBreaks=c(0.4,0.45,0.6,0.8,0.85)
	abline(h=RhosBreaks,lty=2,lwd=2)
	
	if(interactiveMode){line <- readline(interactivePrompt)}
	hist(abs(Correlations[upper.tri(Correlations)]),20,col="red",xlab="abs(Rho)",ylab="Frequency",main="Distribution of absolute values of Rho") 
	
	if(interactiveMode){line <- readline(interactivePrompt)}
	hist(-log(PValues[upper.tri(PValues)],10),20,col="red",xlab="-log10(P-Value)",ylab="Frequency",main="Distribution of P-Values (log scale)") 
	abline(v=-log(0.05,10),lty=2,lwd=2)
	
	if(interactiveMode){line <- readline(interactivePrompt)}
	plot(-log(PValues[upper.tri(PValues)],10)~abs(Correlations[upper.tri(Correlations)]),pch=21,col="black",bg="red",xlab="abs(Rho)",ylab="-log10(P-Value)",main="P-Values as a function of Rhos") 
	abline(h=-log(0.05,10),lty=2,lwd=2)
	RhoMin=abs(Correlations[which.min(abs(PValues-0.05))])
	abline(v=RhoMin,lty=2,lwd=2)
	
	if(interactiveMode){line <- readline("Press return to close the plot...");dev.off()}
}
