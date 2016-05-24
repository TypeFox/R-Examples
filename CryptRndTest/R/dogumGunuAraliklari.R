dogumGunuAraliklari <-function(e,m){
	ee=0
	n=length(e)
	dat=matrix(e[1:(m*floor(n/m))],nrow=floor(n/m),byrow=TRUE)  	
    dat.sorted=unlist(t(apply(dat,1,sort)) )
	dat.diff=t(diff(t(dat.sorted)))
	dat.num.duplicated=apply(dat.diff,1,function(dat.diff){sum(duplicated(dat.diff)==TRUE)}) 

	result=list(ee=dat.num.duplicated,n=n)
	return(result)
}
