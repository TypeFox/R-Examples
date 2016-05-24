logLin=function(data,count.columns,k=10,zero.na=FALSE) {
	nn=data[,count.columns]
	sums=apply(nn,1,function(x){sum(x,na.rm=T)})
	msum=mean(sums)/sums
	nn=nn*msum
	vtr=c()
	for (n in 1:length(nn[,1])) {
		nx=nn[n,]
		nxk=(nx<=k)
		nxK=(nx>k)
		nx[nxK]=log(nx[nxK])
		nx[nxk]=nx[nxk]/k+log(k)-1
		vtr=data.frame(rbind(vtr,nx))
	}
	if (zero.na){
		vtr[nn==0]=NA
	}
	return(vtr)
}
