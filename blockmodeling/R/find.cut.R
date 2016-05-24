"find.cut" <-
function(
	M,	#matrix of a network
	clu,	#partition
	alt.blocks="reg", #alternative block to null block
	cuts="all", #maximumvnumber of evaluations at different cuts
  ... #other parameters to crit.fun
){
  if(cuts=="all"){
    allvals<-sort(unique(M))
  #  allvals<-allvals[allvals>0]
    if(length(allvals)>1000) cat(length(allvals), "evaluations will be made.\n")
    cuts<-allvals
  }

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
  res.IM<-array(NA,dim=c(k2[1],k2[2],length(cuts)))
  res.IM[,,1]<-alt.blocks
  for(i in 1:length(cuts)) res.IM[,,i]<-crit.fun(M=M,clu=clu,blocks=c("null",alt.blocks),cut=cuts[i],approach="bin",...)$IM
  cut<-matrix(NA,nrow=k2[1],ncol=k2[2])
  for(i in 1:k2[1]){
    for(j in 1:k2[2]){
      cut[i,j]<- max(cuts[which(res.IM[i,j,]==alt.blocks)])
    }
  }
  return(cut)
}

