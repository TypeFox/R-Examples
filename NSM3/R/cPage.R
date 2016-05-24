cPage<-function(alpha, k, n, method=NA, n.mc=10000){
  outp<-list()
  outp$stat.name<-"Page L"
  outp$n.mc<-n.mc  
  
  
  if(alpha>1||alpha<0||class(alpha)!="numeric"){
    cat('Error: Check alpha value! \n')
    return(alpha)
  }     

  outp$alpha<-alpha
  outp$n<-n 
  outp$k<-k 
  outp$n.mc<-n.mc
    
  
  ##When the user doesn't give us any indication of which method to use, try to pick one.
  if(is.na(method)){
    if(factorial(outp$k)^outp$n<=10000){
      method<-"Exact"
    }
    if(factorial(outp$k)^outp$n>10000){
      method<-"Monte Carlo"
    }
  }
  #####################################################################
  
  outp$method<-method
  
  L.calc<-function(x){
    return(sum((1:outp$k)*colSums(t(apply(x,1,rank)))))    
  }
  
  possible.ranks<-matrix(rep(1:outp$k,outp$n),ncol=outp$k,byrow=T)
  
  if(outp$method=="Exact"){
    
    possible.perm<-multCh7(possible.ranks)
    exact.dist<-numeric(factorial(outp$k)^outp$n)
    for(i in 1:factorial(outp$k)^outp$n){
      exact.dist[i]<-L.calc(possible.perm[,,i])
    }
    
    
      L.vals<-sort(unique(exact.dist))
      L.probs<-as.numeric(table(exact.dist))/(factorial(outp$k)^outp$n)
      L.dist<-cbind(L.vals,L.probs)
      upper.tails<-cbind(rev(L.dist[,1]),cumsum(rev(L.dist[,2])))
      outp$cutoff.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),1]
      outp$true.alpha.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),2]  
  }

if(outp$method=="Monte Carlo"){
  mc.perm<-matrix(ncol=outp$k,nrow=outp$n)
  mc.stats<-numeric(n.mc)
  for(i in 1:n.mc){
    for(j in 1:n){
      mc.perm[j,]<-sample(possible.ranks[j,])
    }
    mc.stats[i]<-L.calc(mc.perm)
  }
  
  mc.vals<-sort(unique(mc.stats))
  mc.dist<-as.numeric(table(mc.stats))/n.mc
  
  upper.tails<-cbind(rev(mc.vals),cumsum(rev(mc.dist)))
  outp$cutoff.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),1]
  outp$true.alpha.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),2]
  
}

if(outp$method=="Asymptotic"){
  outp$stat.name<-"Page L*"
  outp$cutoff.U<-qnorm(1-alpha)
}

class(outp)<-"NSM3Ch7c"
  return(outp)
}