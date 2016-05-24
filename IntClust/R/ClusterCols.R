ClusterCols <- function(x,Data,nrclusters=NULL,cols=NULL,ColorComps=NULL) {
	
	if(is.null(nrclusters) & is.null(ColorComps)){
		return(x)
	}
	

	else if(!is.null(nrclusters)){
		if(length(cols)<nrclusters){
			stop("Not for every cluster a color is specified")
		}
	}	
	
	if(!is.null(nrclusters)){
		Clustdata=stats::cutree(Data,nrclusters)
		Clustdata=Clustdata[Data$order]

		ordercolors=Clustdata
		order=seq(1,nrclusters)
	
		for (k in 1:length(unique(Clustdata))){
			select=which(Clustdata==unique(Clustdata)[k])
			ordercolors[select]=order[k]
		}
	}
	else{
		cols=rep("black",length(Data$order.lab))
		names(cols)=Data$order.lab
		cols[which(names(cols)%in%ColorComps)]="#EE1289"
		
	}
	
	colfunc=function(x,cols,ColorComps){
		if(is.null(ColorComps)){
		
			indextemp=which(attr(Data$diss,"Labels")==x)
			index1=which(Data$order==indextemp)	
		
			index2=ordercolors[index1]
		
			color=cols[index2]

		}
		else{
			color=cols[which(names(cols)==x)]
		}
		
		return(color)	
	}
	
	if (stats::is.leaf(x)) {
		## fetch label
		label <- attr(x, "label") 
		## set label color to clustercolor
		attr(x, "nodePar") <- list(pch=NA,lab.col=colfunc(label,cols,ColorComps),lab.cex=0.9,font=2)
		attr(x, "edgePar") <- list(lwd=2,col=colfunc(label,cols,ColorComps))
	}
	return(x)
}
