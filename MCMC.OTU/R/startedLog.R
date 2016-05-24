startedLog=function(data,count.columns,logstart=0.1) {
	nn=data[,count.columns]
	sums=apply(nn,1,function(x){sum(x,na.rm=T)})
	msum=mean(sums)/sums
	nn=nn*msum
	logs=log(nn+logstart)
	if(logstart==0) {
		logs[nn==0]=NA
	}
	return(logs)
}
