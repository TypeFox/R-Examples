hyp.group.dist.m.f <-
function(dat,grouping,dist.method="jaccard",distances=FALSE){
	type.c="centroid"
	if(is.numeric(grouping)==TRUE){
		grouping[grouping==0]<-"NO"
		grouping[grouping==1]<-"YES"
			}
		dat.by.host<-t(dat)
		dismatrix<-vegdist(dat.by.host,method=dist.method,binary=TRUE)
		tempdisper<-betadisperF(dismatrix,group=grouping,type=type.c)
		tot.breadth<-sum(tempdisper$distances[which(grouping=="YES")])
		if(distances==FALSE){return(tot.breadth)}else{return(list(tot.breath=tot.breadth,distances.all=tempdisper$distances,distances=tempdisper$distances[which(grouping=="YES")]))}}
