`mat.h` <-
function(y,x,noanalogues, geodist, thresh){
  if(!inherits(y,"dist")){
    if(is.data.frame(y)||!(ncol(y)==nrow(y)&sum(diag(y))==0)){ 
      y<-as.matrix(dist(sqrt(y)))
    }
  }
  y<-as.matrix(y)
  diag(y)<-Inf
  if(inherits(geodist,"dist"))geodist=as.matrix(geodist)

    sapply(1:nrow(y),function(n){
     exneigh<-geodist[n,]>=thresh
     mean(x[exneigh][which(rank(y[n,][exneigh], ties.method="random")<=noanalogues)])
   })
}

