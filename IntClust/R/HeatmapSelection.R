HeatmapSelection<-function(Data,type=c("data","dist","clust","sim"),distmeasure="tanimoto",normalize=FALSE,method="Q",cutoff=NULL,percentile=FALSE,dendrogram=NULL,width=7,height=7){
	
	#create binary similarity heatmap first
	if(type=="data"){
		ClustData<-Cluster(Data=Data,distmeasure=distmeasure,normalize=normalize,method=method,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
		Data=ClustData$DistM
		type="dist"
	}
	
	
	if(type=="clust"){
		Dist=Data$DistM
		if(0<=min(Dist) & max(Dist)<=1){
			SimData=1-Dist
		}
		else{
			NormData=Normalization(Dist,method="Range")
			SimData=1-NormData
		}
		if(is.null(dendrogram)){
			dendrogram=Data
		}
	}
	
	else if(type=="dist"){
		if(0<=min(Data) & max(Data)<=1){
			SimData=1-Data
			if(is.null(dendrogram)){
				dendrogram=Cluster(Data=Data,type="dist",distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
			}
		}
		else{
			NormData=Normalization(Data,method="Range")
			SimData=1-NormData
			if(is.null(dendrogram)){
				dendrogram=Cluster(Data=Data,type="dist",distmeasure="tanimoto",normalize=TRUE,method="Q",clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
			}
		}
		
		
	}
	else if(type=="sim"){
		SimData=Data
		if(0<=min(SimData) & max(SimData)<=1){
			if(is.null(dendrogram)){
				DistData=1-Data
				ClustData=Cluster(Data=DistData,type="dist",distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
			}
		}
		else{
			if(is.null(dendrogram)){
				NormData=Normalization(Dist,method="Range")
				DistData=1-Data
				ClustData=Cluster(Data=DistData,type="dist",distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)		
			}
		}
	}	
	
	
	if(!is.null(cutoff)){
		if(percentile==TRUE){
			cutoff=stats::quantile(SimData[lower.tri(SimData)], cutoff)
		}
		
		SimData_bin <- ifelse(SimData<=cutoff,0,SimData) # Every value higher than the 90ieth percentile is kept, all other are put to zero
	}
	
	else{
		SimData_bin=SimData
	}
	
	
	
	dend <- stats::as.dendrogram(dendrogram$Clust)
	Ind <- stats::order.dendrogram(dend)
	
	SimData_bin=SimData_bin[Ind,Ind]
	
	#Layout<-rbind(4:3, 2:1)
	#lhei <- c(0.4, 4)	
	#lwid <- c(0.4, 4)
	#layout(Layout, widths = lwid, heights = lhei, respect = FALSE)
	grDevices::dev.new(width=width,height=height)
	graphics::par(mar = c(9,7, 7, 9))
	graphics::image(x=1:nrow(SimData_bin),y=1:ncol(SimData_bin),z=t(SimData_bin),col=(grDevices::gray(seq(0.9,0,len=1000))),axes=FALSE,xlab="",ylab="")
	graphics::axis(1, 1:ncol(SimData_bin), labels = colnames(SimData_bin), las = 2, line =0, tick = 0, cex.axis = 0.6)
	graphics::axis(4, 1:nrow(SimData_bin), labels = rownames(SimData_bin), las = 2, line = 0, tick = 0, cex.axis = 0.6)
	
	points=graphics::locator(n=2,type="l")	
	cols=c(floor(points$x[1]),ceiling(points$x[2]))
	rows=c(floor(points$y[1]),ceiling(points$y[2]))
	
	if(cols[1]>cols[2]){
		colseq=seq(cols[2],cols[1],1)
	}
	else{
		colseq=seq(cols[1],cols[2],1)
	}
	
	if(rows[1]>rows[2]){
		rowseq=seq(rows[2],rows[1],1)
	}
	else{
		rowseq=seq(rows[1],rows[2],1)
	}
				
	
	SubsetData=SimData_bin[rowseq,colseq]
	DelRows=rownames(SubsetData)[which(rowSums(SubsetData)==1)]
	DelCols=colnames(SubsetData)[which(colSums(SubsetData)==1)]
	
	if(length(DelRows)!=0 & length(DelCols)!=0){
		Subset=SubsetData[-which(rownames(SubsetData)%in%c(DelRows,DelCols)),-which(colnames(SubsetData)%in%c(DelRows,DelCols))]
	}
	else if(length(DelRows)!=0 & length(DelCols)==0){
		Subset=SubsetData[-which(rownames(SubsetData)%in%c(DelRows)),]
		
	}
	else if(length(DelRows)==0 & length(DelCols)!=0){
		Subset=SubsetData[,-which(colnames(SubsetData)%in%c(DelCols))]
		
	}	
	else if(length(DelRows)==0 & length(DelCols)==0){
		Subset=SubsetData
		
	}
	SelComps=colnames(Subset)
	
	return(SelComps)
}
