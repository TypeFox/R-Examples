### this is the function for selecting the media from multiple optimization results ###

selection<-function(M){
  q<-vector()
  for(i in 1:13){
    p<-vector()
    for(j in 1:length(M)){
      p[j]<-M[[j]][[2]][1,i]
    }
    q<-rbind(q,p,deparse.level=0)
  }
  m<-vector()
  for(i in 1:13){
    m[i]<-as.integer(median(q[i,]))
  }
  return(m)
}