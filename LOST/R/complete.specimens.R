complete.specimens <-
function(dataset,nlandmarks){
  base<-c(1,1)
  included<-base
  excluded<-base
  nspecimen<-nrow(dataset)/nlandmarks
  start<-seq(from=1,to=nrow(dataset),by=nlandmarks)
  nsp<-1:nspecimen
  for (k in 1:nspecimen){
    x<-start[k]
    y<-x+nlandmarks-1
    single<-dataset[x:y,]
    reduced<-na.omit(single)
    rows<-nrow(reduced)
    if (rows==nlandmarks){included<-rbind(included,single)}
    else {excluded<-rbind(excluded,single)}
  }
  end<-nrow(included)
  included<-included[2:end,]
  return (included)
}
