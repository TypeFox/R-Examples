index.G1d<-function(d,cl)
{
  d<-as.matrix(d)
  wgss<-0
  bgss<-0
  for(i in 1:(nrow(d)-1))
  for(j in (i+1):nrow(d)){
    if(cl[i]==cl[j]){
      wgss<-wgss+d[i,j]^2
    }
    else{
      bgss<-bgss+d[i,j]^2
    }
  }
  (bgss/(max(cl)-1))/((sum(d^2)/2-bgss)/(length(cl)-max(cl)))
}

