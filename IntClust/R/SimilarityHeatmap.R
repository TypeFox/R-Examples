SimilarityHeatmap<-function(Data,type=c("data","clust","sim","dist"),distmeasure="tanimoto",normalize=FALSE,method="Q",cutoff=NULL,percentile=FALSE,plottype="new",location=NULL){
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
		
		ClustData=Data
	}
	
	else if(type=="dist"){
		if(0<=min(Data) & max(Data)<=1){
			SimData=1-Data
			ClustData=Cluster(Data=Data,type="dist",distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
			
		}
		else{
			NormData=Normalization(Data,method="Range")
			SimData=1-NormData
			ClustData=Cluster(Data=Data,type="dist",distmeasure="tanimoto",normalize=TRUE,method="Q",clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
			
		}
		
		
	}
	else if(type=="sim"){
		SimData=Data
		if(0<=min(SimData) & max(SimData)<=1){
			DistData=1-Data
			ClustData=Cluster(Data=DistData,type="dist",distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
			
		}
		else{
			NormData=Normalization(Dist,method="Range")
			DistData=1-Data
			ClustData=Cluster(Data=DistData,type="dist",distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)		
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
	
	plottypein(plottype,location)
	gplots::heatmap.2(SimData_bin, 
			Rowv = stats::as.dendrogram(stats::as.hclust(ClustData$Clust)), Colv=stats::as.dendrogram(stats::as.hclust(ClustData$Clust)),trace="none",
			col=(grDevices::gray(seq(0.9,0,len=1000))),
			cexRow=0.6, cexCol=0.6,  
			margins=c(9,9),
			key=FALSE,
			keysize=0.4,
			symkey=FALSE,
			sepwidth=c(0.01,0.01),
			sepcolor="black",
			colsep=c(0,ncol(SimData_bin)),
			rowsep=c(0,nrow(SimData_bin))
	)
	plottypeout(plottype)
	
}
