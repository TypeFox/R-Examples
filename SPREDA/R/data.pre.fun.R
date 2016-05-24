data.pre.fun <-
function(dat){

  ids=unique(dat[,1])
  nn=length(ids)
  res=NULL
  for(i in 1:nn){
    idx=(dat[,1]==ids[i])
    dat.i=dat[idx,]
    dat.i[,2]=dat.i[,2]-dat.i[1,2]+1
    res=rbind(res, dat.i)
  }  
  return(res)
}
