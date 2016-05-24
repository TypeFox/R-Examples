xmat.obj.to.xmat <-
function(dat, dyn.dat, id, time){
  ids1=unique(dat[, id])
  ids2=unique(dyn.dat[, id])
  index.tmp=match(c(id, time), colnames(dyn.dat), 0L)
  res=NULL
  if(sum(!(ids1 %in% ids2))>0){
    print("The id in observation data can not be found in dynamic data")
  }else{
    for(i in 1:length(ids1)){
      idx1=(dat[,id]==ids1[i])
      idx2=(dyn.dat[,id]==ids1[i])
      dat.i=dat[idx1,]
      dyn.dat.i=dyn.dat[idx2,]
      index=findInterval(dat.i[,time], dyn.dat.i[,time])
      tmp=cbind(dat.i, dyn.dat.i[index, -index.tmp])
      res=rbind(res, tmp)
    }
  }
  return(res)
}
