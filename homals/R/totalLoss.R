totalLoss<-function(dframe,x,y,active,rank,level,sets){
  nobj<-nrow(x); ndim<-ncol(x); nset<-length(sets);	loss=0.0
  for (l in 1:nset) {
  	indi<-sets[[l]]; jndi<-indi[which(active[indi])]
    if (length(jndi) == 0) next()
  	ss<-sumSet(dframe,nobj,ndim,y,jndi)
    ii<-which(!is.na(dframe[,jndi[1]]))
	loss<-loss + sum((x[ii,]-ss[ii,])^2)
	}
return(loss)
}
