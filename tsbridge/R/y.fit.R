y.fit <-
function(bug,sims,ysim=NULL,pre.beg=FALSE){
  if(is.null(colnames(sims)))
    stop("columns of sims can not be NULL. names should correspond to parameters")
  if(class(bug)!="tsbugs")
    stop("bug must be a object of class tsbugs")
  
  k<-nrow(sims)
  y<-bug$data$y
  n<-bug$info$n
  beg<-bug$info$args$beg
  temp<-theta.it(bug,sims)
  phi<-temp$phi
  max.phi<-ncol(phi)-1
  ymean<-matrix(NA,k,n)
  ylag<-matrix(0,k,max.phi)
  for(t in beg:n){ 
    ylag[]<-rep(y[(t-1):(t-max.phi)],each=k)
    #missing y use ysim. asked about this on stack exchange
    if(sum(is.na(ylag))>0){
      if(is.null(ysim))
        stop("missing y in data, need some ysim")
      ysimlag<-ysim[(t-1):(t-max.phi)-beg+1]
      ylag[is.na(ylag)]<-ysimlag[is.na(ylag)]
    }
    ymean[,t]<-phi[,1]+rowSums(phi[,-1]*(ylag-phi[,1]))
  }
  
  temp<-ymean
  if(pre.beg==FALSE)  temp<-ymean[,-(1:(beg-1))]
  return(temp)
}
