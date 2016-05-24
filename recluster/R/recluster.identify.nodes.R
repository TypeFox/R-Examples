recluster.identify.nodes<-function (mat,low=TRUE) {
	mat<-as.matrix(mat)
	mat2<-cbind(mat, c(1:nrow(mat)))
	rownames(mat2)<-mat2[,ncol(mat2)]
	mat3<-mat2
	mat2<-mat2[complete.cases(mat2),]
	mat2<-mat2[-1,]
	out<-NULL
	pamrun<-NULL
	val<-NULL
	sil<-NULL
	diffmean<-NULL
	for (k in 1 :(ncol(mat2)-1)){
		pamrun<-pam(mat2[,k],2,diss=F, metric = "euclidean")
		sil[k]<-pamrun$silinfo$avg.width
		type<-pamrun$clustering
		mean<-aggregate(mat2[,k], by=list(type), FUN=mean)
		mean1<-aggregate(mat2[,k], by=list(type), FUN=mean)[,2]
		diffmean[k]<-abs(mean1[1]-mean1[2])
		if(low=="TRUE"){mult<-ncol(mat2)-k}
		else {mult<-1}
		val[k]<-sil[k]*diffmean[k]*mult	
		}
	best<-which.max(val)
	best_pam<-pam(mat2[,best],2,diss=F, metric = "euclidean")
	type<-best_pam$clustering
	mean<-aggregate(mat2[,best], by=list(type), FUN=mean)
	if(mean[1,2] < mean[2,2]){type<-abs(type-3)}	
	x<-function (x)
		{a<-aggregate(x, by=list(type), FUN=mean)[,2]}
	mean2<-apply(mat2[,1:(ncol(mat2)-1)],2,x)	
	mat2<-cbind(mat2,type)
	matrix<-merge(mat2[,ncol(mat2)], mat3[,ncol(mat3)], by = "row.names", all = TRUE)
	matrix<-matrix[order(as.numeric(matrix[,1])),]
	out$mean<-mean2
	out$nodes<-matrix[,2]
	out$scale<-best
	tr<-t(mat2)
	plot(tr[1:(nrow(tr)-2),1], type="l", ylim=c(0,100),xlab="Step",ylab="Support",col=type[1])
	for (i in 2:ncol(tr)){
		lines(tr[1:(nrow(tr)-2),i], type="l",col=type[i])
		}
	a<-c(1:ncol(out$mean))
	points(a,t(out$mean)[,1],col=1,pch=18)
	points(a,t(out$mean)[,2],col=2,pch=18)
	return(out)
}
