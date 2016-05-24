Ptensor<-function(x, k){

  rankx<-dim(x)[2]
  if(is.null(rankx)==FALSE){
    x<-scale(x, scale=FALSE)
    A<-to.tensor(1:(rankx^k), rep(rankx,k))
    comb.pos<-expand.grid(lapply(1:k, function(x){1:rankx}))
    for(i in 1:length(comb.pos[,1])){
      A[i]<-sum(apply(x[,unlist(comb.pos[i,])], 1, prod))/dim(x)[1]
    }
    return(A)
  }else{
    if(k==1){
      return(mean(x))
    }else{
      return(mean((x-mean(x))^k))
    }
  }
}

