sig.number<-function(fdr.table,FDR=0.05,qt=-1){
	if (qt==0.25)
	k<-2
	if  (qt==0.5)
	k<-3
	if  (qt==0.75)
	k<-4
	if   (qt==0.9)
	k<-6
	if   (qt==-1)
	k<-6
	temp<-rep(FDR,times=dim(fdr.table)[1])
	r<-which(fdr.table[,k]<=temp)
	if (length(r)>0){
		ix<-r[length(r)]
	}
	if (length(r)==0){
		ix<-1
	}
	return (fdr.table[ix,1])
}