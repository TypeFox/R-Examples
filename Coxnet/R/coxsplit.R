

########################
#####  Split data  #####
########################

coxsplit=function(y, nfolds){
  N=nrow(y)
  foldid=sample(rep(seq(nfolds), length=N))
  return(foldid)
}

coxsplity=function(y, nfolds){
  N=nrow(y)
  tem=data.frame(y, i=seq(N), foldid=0)
  tem=tem[order(y[, "time"], y[, "status"]), ]
  n1=sum(y[, "status"]);n2=N-n1
  
  tem$foldid[tem[, "status"]==1]=sample(rep(seq(nfolds), length=n1))
  tem$foldid[tem[, "status"]==0]=sample(rep(seq(nfolds), length=n2))
  
  foldid=tem$foldid[order(tem$i)]
  return(foldid)
}

coxsplitw=function(w, nfolds){
  N=length(w)
  tem=data.frame(w=w,i=rep(1:nfolds,time=ceiling(N/nfolds))[1:N])
  tem=tem[order(tem$w),]
  tem=data.frame(tem,foldid=rep(1:nfolds,time=ceiling(N/nfolds))[1:N])
  foldid=tem$foldid[order(tem$i)]
  return(foldid)
}


