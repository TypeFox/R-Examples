coeplot <-
function(coe.lm,obj.index,ylim=c(-1,1)){

	col.vec<-c("red","yellowgreen","green","blue")

	index.mat<-obj.index$mat
	index.flat<-obj.index$flat

	coe.mat<-index.mat

	for(i in 1:(nrow(coe.lm)-1)){
		pos_row<-which(rownames(coe.mat)==index.flat[i,2])
		pos_col<-which(colnames(coe.mat)==index.flat[i,3])
		coe.mat[pos_row,pos_col]<-coe.lm[i+1,1]
	}

	coe.mat[is.na(coe.mat)]<-0

	axi.y<-as.numeric(formatC(seq(ylim[1],ylim[2],length=7),digits=2))

	plot(coe.mat[1,],col=col.vec[1],type="l",lwd=2,ann=FALSE,axes=FALSE,ylim=ylim)
	axis(1, at=1:length(coe.mat[1,]),labels=colnames(coe.mat), col.axis="black")
	axis(2, at=axi.y,labels=axi.y, col.axis="black")
	for(j in 2:nrow(coe.mat)){
		lines(coe.mat[j,],col=col.vec[j],lwd=2)
	}
	title(main=NULL,xlab="Pos",ylab="Coefficient")
	legend(1, ylim[2], rownames(coe.mat), col = col.vec, lty=1, lwd=2,cex = 1)
}
