missing.specimens <-
function(dataset,nspremove,nldremove,nlandmarks){
  
  remove.coor.random<-function(single,j,nlandmarks){
    nld<-1:nlandmarks
    r<-sample(nld,j)
    for (m in 1:j){
      d<-r[m]
      single[d,]<-c(NA,NA)
    }
    return(single)
  }
  
  nspecimen<-nrow(dataset)/nlandmarks
  start<-seq(from=1,to=nrow(dataset),by=nlandmarks)
  nsp<-1:nspecimen
  removes<-sample(nsp,nspremove)
  for (k in 1:nspremove){
    specimen<-removes[k]
    x<-start[specimen]
    y<-x+nlandmarks-1
    single<-dataset[x:y,]
    j<-sample(nldremove,1)
    newsingle<-remove.coor.random(single,j,nlandmarks)
    dataset[x:y,]<-newsingle}
  return(dataset)
}
