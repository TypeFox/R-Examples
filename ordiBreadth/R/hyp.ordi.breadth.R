hyp.ordi.breadth <-
function(dat,grouping,dist.method="jaccard",distance=FALSE){
	type="centroid"
	if(is.numeric(grouping)==TRUE){
		grouping[grouping==0]<-FALSE
		grouping[grouping==1]<-TRUE
			}
		dat.by.host<-t(dat)
		dismatrix<-vegdist(dat.by.host,method=dist.method,binary=TRUE)
		tempdisper<-betadisperF(dismatrix,group=grouping,type=type)
		tot.breadth<-sum(tempdisper$distances[which(grouping==TRUE)])
		if(dim(tempdisper$centroids)[1]==1){rownames(tempdisper$centroids)=TRUE}
		if(distance==FALSE){return(tot.breadth)}else{return(list(tot.breadth=tot.breadth,distances=tempdisper$distances[which(grouping=="TRUE")],centroid=tempdisper$centroids[which(rownames(tempdisper$centroids)==TRUE),]))}
		
	}
