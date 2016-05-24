firstValueRow <-
function(x) {
	cumSumMat<-matrix(NA, nrow=dim(x)[1], ncol=dim(x)[2])
	for(i in 1:(dim(x)[1])) {
		cumSumMat[i,]<-cumsum(x[i,])
	}
	cumSumMat2<-cbind(matrix(0, nrow=dim(x)[1], ncol=1), cumSumMat[,-(dim(x)[2])])
	ResultMat<-matrix(NA, nrow=dim(x)[1], ncol=dim(x)[2])
	for(i in 1:dim(x)[2]) {
		ResultMat[,i]<-ifelse(cumSumMat2[,i]>0, 0, x[,i])
	}
	return(ResultMat)	
}
