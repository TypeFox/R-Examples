ComparePlot<-function(List,nrclusters=NULL,cols=NULL,fusionsLog=FALSE,WeightClust=FALSE,names=NULL,margins=c(8.1,3.1,3.1,4.1),plottype="new",location=NULL,...){
	
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	
	for(i in 1:length(List)){
		if(attributes(List[[i]])$method == "Weighted" & WeightClust==TRUE){
			T=List[[i]]$Clust
			attr(T,"method")="Single Clustering"
			List[[i]]=T
		}
	}
	
	MatrixColors=ReorderToReference(List,nrclusters,fusionsLog,WeightClust,names)
	
	Names=ColorsNames(MatrixColors,cols)
	
	nobs=dim(MatrixColors)[2]
	nmethods=dim(MatrixColors)[1]
	
	if(is.null(names)){
		for(j in 1:nmethods){
			names[j]=paste("Method",j,sep=" ")	
		}
	}
	
	#similar=round(SimilarityMeasure(MatrixColors),2)
	plottypein(plottype,location)
	graphics::par(mar=margins)
	plotrix::color2D.matplot(MatrixColors,cellcolors=Names,show.values=FALSE,axes=FALSE,xlab="",ylab="",...)
	graphics::axis(1,at=seq(0.5,(nobs-0.5)),labels=colnames(MatrixColors),las=2,cex.axis=0.70)
	graphics::axis(2,at=seq(0.5,(nmethods-0.5)),labels=rev(names),cex.axis=0.70,las=2)
	#axis(4,at=seq(0.5,(nmethods-0.5)),labels=rev(similar),cex.axis=0.65,las=2)
	plottypeout(plottype)
}
