.kfwer.distrib <- function(permT, k=1){
  #permT Bxm matrix of test statistics (NB: larger is better)
  permT=t(apply(permT,1,sort,decreasing = TRUE))
  if(k>ncol(permT)) {
    k=ncol(permT)
    warning(paste("k=",k," is greater than the total number of hypotheses. k is now set to ",ncol(permT),sep=""))    
  }
  wk=permT[,k,drop=FALSE]
}

.p.adjust.kfwer <- function(permT,k=1){
  wk=.kfwer.distrib(permT=permT,k=k)
  adj.ps=sapply(permT[1,],function(x) (sum(x<=wk)+1)/(length(wk)+1))
  adj.ps
}