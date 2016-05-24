`mergeBlockup` <-
function(blocklist,blockvalues,x,w,i,block){
  n<-length(blockvalues); nn<-1:n; ii<-which(i+1!=nn)
  blocklist[i,]<-c(blocklist[i,1],blocklist[i+1,2])
  indi<-blocklist[i,1]:blocklist[i+1,2]
  blockvalues[i]<-block(x[indi],w[indi])
  blocklist<-blocklist[ii,]
  if (length(ii) == 1) dim(blocklist)<-c(1,2)
  blockvalues<-blockvalues[ii]
  list(v=blockvalues,l=blocklist)
}

