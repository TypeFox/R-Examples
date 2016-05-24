match.dat.fun <-
function(dat1, dat2){
  ids=as.matrix(unique(dat1[, 1]))
  n.len=length(ids)
  res=NULL
  for(i in 1:n.len){
    idx=(dat2[,1]==ids[i])
    dat.i=dat2[idx, ]
    res=rbind(res, dat.i)
  }
  return(res)
}
