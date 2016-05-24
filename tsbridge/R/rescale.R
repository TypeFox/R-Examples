rescale <-
function(bug, sims, to.real=NULL){
  if(is.null(colnames(sims)))
    stop("columns of sims can not be NULL. names should correspond to parameters")
  if(class(bug)!="tsbugs")
    stop("bug must be a object of class tsbugs")
  
  theta<-colnames(sims)
  stoc<-NULL
  dd<-nodes(bug, part="prior")
  dd<-subset(dd, stoc==1)
  for(i in 1:length(theta)){
    if(dd$stoc[i]!=1)
      stop("use only give variables in sims that are from stochastic nodes, not transformed nodes")
    if(dd$dist[i]=="dgamma"){
      if(to.real==TRUE){
        sims[,i]<-log(sims[,i])
      }
      if(to.real==FALSE){
        sims[,i]<-exp(sims[,i])
      }
    }
    f1<-function(x,p1=min(x),p2=max(x)){
      y<-(x-p1)/(p2-p1)
      y
    }
    f2<-function(x,p1=min(x),p2=max(x)){
      y<-x*(p2-p1)+p1
      y
    }
    if(dd$dist[i]=="dunif"){
      if(to.real==TRUE){
        sims[,i]<-qlogis(f1(sims[,i],dd$param1[i]+0.0001,dd$param2[i]+0.0001))
      }
      if(to.real==FALSE){
        sims[,i]<-f2(plogis(sims[,i]),dd$param1[i]+0.0001,dd$param2[i]+0.0001)
      }
    }
    if(dd$dist[i]=="dbeta"){
      if(to.real==TRUE){
        sims[,i]<-qlogis(f1(sims[,i],p1=0,p2=1))
      }
      if(to.real==FALSE){
        sims[,i]<-f2(plogis(sims[,i]),p1=0,p2=1)
      }
    }
  }
  sims
}
