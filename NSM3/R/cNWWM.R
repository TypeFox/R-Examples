cNWWM<-function(alpha, k, n, method=NA, n.mc=10000){
  outp<-list()
  outp$stat.name<-"Nemenyi, Wilcoxon-Wilcox, Miller R*"
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
  
  R.star.calc<-function(two.dim.mat,u){
    row.ranks<-t(apply(two.dim.mat,1,rank))
    return(colSums(row.ranks)[u]-colSums(row.ranks)[1])
  }
  
  R.star.all<-function(two.dim.mat){
    row.ranks<-t(apply(two.dim.mat,1,rank))
    return(max(colSums(row.ranks)[-1])-colSums(row.ranks)[1])      
  }
  
  outp$num.comp<-num.comp<-outp$k-1 
  
  possible.ranks<-matrix(rep(1:outp$k,outp$n),ncol=outp$k,byrow=T)
  
  if(outp$method=="Exact"){
    possible.perm<-multCh7(possible.ranks)
    
    exact.dist<-apply(possible.perm,3,R.star.all)
    
    R.star.vals<-sort(unique(exact.dist))
    R.star.probs<-as.numeric(table(exact.dist))/(factorial(outp$k)^outp$n)
    R.star.dist<-cbind(R.star.vals,R.star.probs)
    upper.tails<-cbind(rev(R.star.dist[,1]),cumsum(rev(R.star.dist[,2])))
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
      mc.stats[i]<-R.star.all(mc.perm)
    }
    
    mc.vals<-sort(unique(mc.stats))
    mc.dist<-as.numeric(table(mc.stats))/n.mc
    
    upper.tails<-cbind(rev(mc.vals),cumsum(rev(mc.dist)))
    outp$cutoff.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),1]
    outp$true.alpha.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),2]
  }  
  
  if(outp$method=="Asymptotic"){
    outp$cutoff.U<-cMaxCorrNor(alpha,outp$k-1,0.5)*(outp$n*outp$k*(outp$k+1)/6)^(1/2)
  }
  class(outp)<-"NSM3Ch7c"
  outp
}
