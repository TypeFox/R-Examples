xmat.to.cumsum <-
function(dat){
  ids=unique(dat[, 1])
  n=length(ids)
  res=NULL

  for(i in 1:n){
    idx=(dat[,1]==ids[i])
    dat.i=dat[idx,]
    time.1=c(0, dat.i[-sum(idx),2])
    time=dat.i[,2]-time.1
    fun.tmp=function(x){
      x*time
    }
    dat.tmp=apply(dat.i[,-c(1,2)], 2, function(x) cumsum(x*time))
    res=rbind(res, dat.tmp)    
  }
  res=as.data.frame( res)
  return(res)
}
