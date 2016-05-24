`minrc` <-
function(x){
	rx<-max(x,na.rm=TRUE)
	n.lev<-numeric(nrow(x))
	for(i in 1:rx){
		tmp<-x==i
		tmp2<-rowSums(tmp,na.rm=TRUE)>0
		n.lev<-n.lev+tmp2
	}
	uni.lev<-sort(unique(n.lev),decreasing=TRUE)
	if(uni.lev[1]!=rx)
		stop("No variable shows all possible levels.",call.=FALSE)
	if(length(uni.lev)==1)
		return(uni.lev)
	mat.rc<-matrix(rx,nrow(x),nrow(x))
	for(i in uni.lev[-1]){
		ids<-which(n.lev==i)
		mat.rc[ids,]<-i
		mat.rc[,ids]<-i
	}
	mat.rc
}

