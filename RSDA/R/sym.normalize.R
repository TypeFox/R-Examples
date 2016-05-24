sym.normalize <-
function(sym.data) {
  nn<-sym.data$N
  mm<-sym.data$M
  pos<-1
  for(j in 1:mm) {
    sdc<-sym.sd(sym.var(sym.data,j))
    mc<-sym.mean(sym.var(sym.data,j))
    for(i in 1:nn) {
      sym.data$meta[i,sym.data$sym.var.starts[j]]<-(sym.data$meta[i,sym.data$sym.var.starts[j]]-mc)/sdc
      sym.data$meta[i,sym.data$sym.var.starts[j]+1]<-(sym.data$meta[i,sym.data$sym.var.starts[j]+1]-mc)/sdc
      sym.data$data[i,pos]<-(sym.data$data[i,pos]-mc)/sdc
      sym.data$data[i,pos+1]<-(sym.data$data[i,pos+1]-mc)/sdc
    }
    pos<-pos+2
  }
  return(sym.data)
}
