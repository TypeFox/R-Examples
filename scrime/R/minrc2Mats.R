`minrc2Mats` <-
function(catX1,catX2,n.cat){
	uni.lev1<-unique(catX1)	
	uni.lev2<-unique(catX2)
	uni.lev<-unique(c(uni.lev1,uni.lev2))
	uni.lev<-sort(uni.lev,decreasing=TRUE)
	if(uni.lev[1]!=n.cat)
		stop("No variable shows all possible levels.",call.=FALSE)
	if(length(uni.lev)==1)
		return(uni.lev)
	mat.rc<-matrix(n.cat,length(catX1),length(catX2))
	for(i in uni.lev[-1]){
		mat.rc[catX1==i,]<-i
		mat.rc[,catX2==i]<-i
	}
	mat.rc
}

