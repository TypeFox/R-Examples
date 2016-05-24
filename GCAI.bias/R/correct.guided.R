correct.guided <-
function(coe.train,obj.test){
	count.n<-obj.test$counts[,2]
	x.ti<-obj.test$index
	coe.t<-coe.train[-1]
	betax<-as.matrix(x.ti)%*%matrix(coe.t,ncol=1)
	count.m<-2^(betax+log2(count.n))

	count.m<-as.integer(count.m+0.5)
	count.m[which(count.m<1)]<-1

	mat.test<-cbind(obj.test$counts,count.m,betax)
	
	colnames(mat.test)[ncol(mat.test)-1]<-paste(colnames(mat.test)[ncol(mat.test)-2],"corrected",sep="_")
	colnames(mat.test)[ncol(mat.test)]<-"betax"

	return(mat.test)	
}
