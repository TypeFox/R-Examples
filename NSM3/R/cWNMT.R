cWNMT<-function(alpha, k, n, method=NA, n.mc=10000){
  outp<-list()
  outp$stat.name<-"Wilcoxon, Nemenyi, McDonald-Thompson R"
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
    if(outp$k*factorial(outp$k)^outp$n<=10000){
      method<-"Exact"
    }
    if(outp$k*factorial(outp$k)^outp$n>10000){
      method<-"Monte Carlo"
    }
  }
  #####################################################################
  
  outp$method<-method  
  
  R.calc<-function(two.dim.mat,u,v){
    row.ranks<-t(apply(two.dim.mat,1,rank))
    return(abs(colSums(row.ranks)[u]-colSums(row.ranks)[v]))
  }
  
  R.all<-function(two.dim.mat){
    row.ranks<-t(apply(two.dim.mat,1,rank))
    return(max(colSums(row.ranks))-min(colSums(row.ranks)))      
  }
  
  outp$num.comp<-num.comp<-outp$k*(outp$k-1)/2  
  
  possible.ranks<-matrix(rep(1:outp$k,outp$n),ncol=outp$k,byrow=T)
  
  if(outp$method=="Exact"){
    possible.perm<-multCh7(possible.ranks)
    
    exact.dist<-numeric(factorial(outp$k)^outp$n)
    for(i in 1:factorial(outp$k)^outp$n){
      exact.dist[i]<-R.all(possible.perm[,,i])
    }
      
    R.vals<-sort(unique(exact.dist))
    R.probs<-as.numeric(table(exact.dist))/(factorial(outp$k)^outp$n)
    R.dist<-cbind(R.vals,R.probs)
    upper.tails<-cbind(rev(R.dist[,1]),cumsum(rev(R.dist[,2])))
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
      mc.stats[i]<-R.all(mc.perm)
    }
    
    mc.vals<-sort(unique(mc.stats))
    mc.dist<-as.numeric(table(mc.stats))/n.mc
    
    upper.tails<-cbind(rev(mc.vals),cumsum(rev(mc.dist)))
    outp$cutoff.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),1]
    outp$true.alpha.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),2]
  }  
  
  if(outp$method=="Asymptotic"){
    outp$cutoff.U<-cRangeNor(alpha,outp$k)*(outp$n*outp$k*(outp$k+1)/12)^(1/2)
  }
  class(outp)<-"NSM3Ch7c"
  outp
}
