cDurSkiMa<-function(alpha,obs.mat, method=NA, n.mc=10000){
  ##Figure out k,n,s,p,lambda based on obs.mat;
  outp<-list()
  outp$stat.name<-"Durbin, Skillings-Mack D"
  outp$n.mc<-n.mc
  
  outp$n<-n<-nrow(obs.mat) 
  outp$k<-k<-ncol(obs.mat) 
  
  outp$pp<-unique(colSums(obs.mat))
  
  if(length(outp$pp)>1){
    stop("Check that p is the same for each treatment.")
  }
  
  outp$ss<-unique(rowSums(obs.mat))
  
  if(length(outp$ss)>1){
    stop("Check that s is the same for each treatment.")
  }
  
  outp$lambda<-outp$pp*(outp$ss-1)/(outp$k-1)
  
  if(alpha>1||alpha<0||class(alpha)!="numeric"){
    cat('Error: Check alpha value! \n')
    return(alpha)
  }
  
  outp$alpha<-alpha
  
  ##When the user doesn't give us any indication of which method to use, try to pick one.
  if(is.na(method)){
    if(outp$ss*factorial(outp$ss)^outp$pp<=10000){
      method<-"Exact"
    }
    if(outp$ss*factorial(outp$ss)^outp$pp>10000){
      method<-"Monte Carlo"
    }
  }
  #####################################################################
  
  outp$method<-method  
  
  outp$obs.mat<-obs.mat
  
  DSK.stat<-function(obs.data){
    tmp.mat<-outp$obs.mat
    for(i in 1:outp$n){
      tmp.mat[i,tmp.mat[i,]!=0]<-obs.data[i,]
    }
    Rj<-apply(tmp.mat,2,function(x) sum(x[!is.na(x)]))
    D.stat<-12/(outp$lambda*outp$k*(outp$ss+1))*sum((Rj-outp$pp*(outp$ss+1)/2)^2)
    return(D.stat)
  }
  
  possible.ranks<-matrix(rep(1:outp$ss,outp$n),ncol=outp$ss,byrow=T)
  
  if(outp$method=="Exact"){
    possible.perm<-multCh7(possible.ranks)
    exact.dist<-apply(possible.perm,3,DSK.stat)
    
    D.vals<-sort(unique(exact.dist))
    D.probs<-as.numeric(table(exact.dist))/(factorial(outp$ss)^outp$n)
    D.dist<-cbind(D.vals,D.probs)
    upper.tails<-cbind(rev(D.dist[,1]),cumsum(rev(D.dist[,2])))
    outp$cutoff.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),1]
    outp$true.alpha.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),2]  
  }
  
  if(outp$method=="Monte Carlo"){
    mc.perm<-matrix(ncol=outp$ss,nrow=outp$n)
    mc.stats<-numeric(n.mc)
    for(i in 1:n.mc){
      for(j in 1:n){
        mc.perm[j,]<-sample(possible.ranks[j,])
      }
      mc.stats[i]<-DSK.stat(mc.perm)
    }
    
    mc.vals<-sort(unique(mc.stats))
    mc.dist<-as.numeric(table(mc.stats))/n.mc
    
    upper.tails<-cbind(rev(mc.vals),cumsum(rev(mc.dist)))
    outp$cutoff.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),1]
    outp$true.alpha.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),2]
  }  
  
  if(outp$method=="Asymptotic"){
    outp$p.val<-qchisq(1-alpha,outp$k-1)    
  }
  class(outp)<-"NSM3Ch7c"
  outp
}
  
  
  