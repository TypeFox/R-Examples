ClusterPlot<-function(Data1,Data2=NULL,nrclusters=NULL,cols=NULL,plottype="new",location=NULL,ColorComps=NULL,...){
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
	
	x=Data1$Clust
	if(is.null(Data2)){
		Data=Data1$Clust
	}
	else{
		Data=Data2$Clust
	}
	
	d_temp<- stats::dendrapply(stats::as.dendrogram(x,hang=0.02),ClusterCols,Data,nrclusters,cols,ColorComps)
	plottypein(plottype,location)
	graphics::plot(d_temp,nodePar=list(pch=NA),edgePar=list(lwd=2),ylab="Height",font.axis=2,font.lab=2,font=2,...)
	graphics::axis(side = 2, lwd = 2)	
	plottypeout(plottype)
}
