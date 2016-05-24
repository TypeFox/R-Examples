h.fit <-
function(bug,sims,pre.beg=FALSE){
  if(is.null(colnames(sims)))
    stop("columns of sims can not be NULL. names should correspond to parameters")
  if(class(bug)!="tsbugs")
    stop("bug must be a object of class tsbugs")
  
  k<-nrow(sims)
  y<-bug$data$y
  n<-bug$info$n
  beg<-bug$info$args$beg
  theta<-theta.it(bug,sims)
  if(bug$info$variance=="SV"){
    psi<-theta$psi
    psi.star<-theta$psi.star
    psi[,1][psi[,1]==0]<--log(psi.star[,1])
    psi[,-1][psi[,-1]==0]<-2*psi.star[,-1]-1
    max.psi<-ncol(psi)-1
    sv.beg<-beg+bug$info$args$sv.order
    #sv.beg<-bug$info$args$sv.beg
    #create h, equivlent to y (data in y.mean.fn)
    h<-matrix(NA, k, n)
    hn<-paste0("h.",1:n)
    colnames(h)<-hn
    hn2<-colnames(theta$h)
    h[,which(hn %in% hn2)]<-theta$h[,hn2[hn2 %in% hn]]
    hmean<-matrix(NA, k, n)
    for(t in beg:(sv.beg-1)){
      hmean[,t]<-psi[,1]
    }
    hlag<-matrix(0,k,max.psi)
    for(t in sv.beg:n){
      hlag[]<-h[,(t-1):(t-max.psi)]
      hmean[,t]<-psi[,1]+rowSums(psi[,-1]*(hlag-psi[,1]))
    }
  }
  if(bug$info$variance=="RV"){
    isig02<-theta$isig02
    beta<-cbind(matrix(NA,k,beg),theta$beta)
    delta<-cbind(matrix(NA,k,beg),theta$delta)
    lsig<-matrix(NA, k, n)
    lsig[,beg]<--0.5*log(isig02)
    for(t in (beg+1):n){
      lsig[,t]<-lsig[,t-1]+(delta[,t]*beta[,t])
    }
    hmean<-2*lsig
  }
  temp<-hmean
  if(pre.beg==FALSE)  temp<-hmean[,-(1:(beg-1))]
  return(temp)
}
