"find.m2" <-
function(
	M,	#matrix of a network
	clu,	#partition
	alt.blocks="reg", #alternative block to null block
	neval=100, #number of evaluations at different ms
	half = TRUE,	# should the returned value of m be one half of the value where the incosnistencies are the same, otherwise, the m is restricted to max(M)
	ms=NULL,	#the values of m where the function should be evaluated
	... #other parameters to crit.fun
){
  if(is.null(ms)){
  	ms<-seq(from=min(M), to=max(M)*(1+half), length.out=neval)
  } else neval<-length(ms)

  if(is.list(clu)){
		k<-sapply(clu,function(x)length(unique(x)))
		clu<-lapply(clu,function(x)as.integer(factor(x)))
		if(length(k)>2) {
			for(i in 2:length(clu)){
				clu[[i]]<-clu[[i]] + max(clu[[i-1]])
  			}
  		k2<-max(clu[[length(clu)]])
  		} else k2<-k
	} else {
		k<-length(unique(clu))
		clu<-as.integer(factor(clu))
		k2<-c(k,k)
	}
  res.IM<-array(NA,dim=c(k2[1],k2[2],length(ms)))
  for(i in 1:neval) res.IM[,,i]<-crit.fun(M=M,clu=clu,blocks=c("null",alt.blocks),m=ms[i],approach="val",...)$IM
  m<-matrix(NA,nrow=k2[1],ncol=k2[2])
  for(i in 1:k2[1]){
    for(j in 1:k2[2]){
      m[i,j]<- max(ms[which(res.IM[i,j,]==alt.blocks)])
    }
  }
  m[m== -Inf]<-0
  if(half) m<-m/2
  return(m)
}
  

