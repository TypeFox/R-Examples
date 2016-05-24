posplot <-
function(test.corrected,pos.test,n.lim=1000, fit.cut.lr=50){

	mat.p<-cbind(pos.test,test.corrected)

	mat.p<-mat.p[order(mat.p[,1], mat.p[,2]),]

	colnames(mat.p)<-gsub("_.*Pos","",colnames(mat.p))

	if(nrow(mat.p)<1000){
		n.lim<-nrow(mat.p)
	}

	par(mfrow=c(4,1))
	par(mar=c(2,4,2,2))
	for(i in 1:3){
		barplot(mat.p[1:n.lim,i+4],ylab=colnames(mat.p)[i+4],ylim=c(0,max(mat.p[1:n.lim,i+4])))
	}

	mat.p<-mat.p[which(rowSums(mat.p[,5:6])>=fit.cut.lr),]

	lr.before<-log2((mat.p[,5]+0.5)/(mat.p[,6]+0.5))
	lr.after<-log2((mat.p[,5]+0.5)/(mat.p[,7]+0.5))

	y.lim=range(c(lr.before, lr.after, 2^mat.p[,ncol(mat.p)]))
	plot(lr.before,ylab="Log Ratio",type="l",ylim=y.lim, lwd=1.5)
	points(2^mat.p[,ncol(mat.p)],type="l",col="red")
	points(lr.after,type="l",col="blue")
}
