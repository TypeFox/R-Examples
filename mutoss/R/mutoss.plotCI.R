############################################
# mutoss.plotCI

mutoss.plotCI<-function(mat){
	diff<-max(mat[,1])-min(mat[,1])
	if(any(is.na(mat[,2])))	mat[,2]<-min(mat[,1])-diff*2
	if(any(is.na(mat[,3])))	mat[,3]<-max(mat[,1])+diff*2
	k<-nrow(mat)
	plotCI(1:k,mat[,1],abs(mat[,3]-mat[,1]),abs(mat[,2]-mat[,1]),lwd=2,col="red",scol="blue",
			main="CI plot",xaxt="n",xlab="Parameters",
			ylab="Values")
	axis(1, at=c(1:k), labels=rownames(mat)) 
}

