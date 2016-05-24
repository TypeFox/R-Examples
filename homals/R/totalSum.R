totalSum<-function(dframe,x,y,active,rank,level,sets){
  nobj<-dim(x)[1]; ndim<-dim(x)[2]; nset<-length(sets);	stot<-array(0.0,dim(x))
  for (l in 1:nset) {
  	indi<-sets[[l]]; jndi<-indi[which(active[indi])]
    if (length(jndi) == 0) next()
  	ss<-sumSet(dframe,nobj,ndim,y,jndi)
    ii<-which(!is.na(dframe[,jndi[1]]))
	stot[ii,]<-stot[ii,] + ss[ii,]
	}
return(stot)
}
