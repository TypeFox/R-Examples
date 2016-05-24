cSkilMack<-function(alpha,obs.mat, method=NA, n.mc=10000){
  outp<-list()
  outp$stat.name<-"Skillings-Mack SM"
  outp$n.mc<-n.mc
  
  outp$obs.mat<-obs.mat
  outp$n<-n<-nrow(obs.mat) 
  outp$k<-k<-ncol(obs.mat) 
  
  if(alpha>1||alpha<0||class(alpha)!="numeric"){
    cat('Error: Check alpha value! \n')
    return(alpha)
  }
  
  outp$alpha<-alpha

  outp$ss<-s<-rowSums(outp$obs.mat)

  ##When the user doesn't give us any indication of which method to use, try to pick one.
  if(is.na(method)){
    if(prod(factorial(outp$ss))<=10000){
      method<-"Exact"
    }
    if(prod(factorial(outp$ss))>10000){
      method<-"Monte Carlo"
    }
  }
  #####################################################################
  
  outp$method<-method
  
  lambda.mat<-matrix(nrow=k,ncol=k)
  for(i in 1:k){
    for(j in (1:k)){
      lambda.mat[i,j]<-sum(outp$obs.mat[,i]*outp$obs.mat[,j])
    }
  }
  
  sigma.mat<-(-1)*lambda.mat[1:(k-1),1:(k-1)]
  for(i in 1:(k-1)){
    diag(sigma.mat)[i]<-sum(lambda.mat[i,-i])
  }
  
  ##Uses MASS##
  sigma0.inv<-ginv(sigma.mat)
  #############

  missing.obs<-function(rank.data){
    si<-sum(!is.na(rank.data))
    rank.data[is.na(rank.data)]<-(si+1)/2
    return(sqrt(12/(si+1))*(rank.data-(si+1)/2))
  }
  
  SM.stat<-function(obs.data){    
    ranks<-t(apply(obs.data,1,rank,na.last="keep"))
    ranks<-t(apply(ranks,1,missing.obs))
    Aj<-apply(ranks,2,sum)[1:(k-1)]
    SM.stat<-t(Aj)%*%sigma0.inv%*%t(t(Aj))
    return(as.numeric(SM.stat))
  }
  
  
  outp$obs.mat[outp$obs.mat==0]<-NA
    
  possible.ranks<-matrix(ncol=outp$k,nrow=outp$n)
  for(i in 1:outp$n){
    possible.ranks[i,!is.na(outp$obs.mat[i,])]<-1:sum(!is.na(outp$obs.mat[i,]))
  }
  
  
  if(outp$method=="Exact"){
    possible.perm<-multCh7SM(possible.ranks)
    exact.dist<-apply(possible.perm,3,SM.stat)
    SM.vals<-sort(unique(round(exact.dist,5)))
    SM.probs<-as.numeric(table(exact.dist))/(prod(factorial(rowSums(!is.na(outp$obs.mat)))))
    SM.dist<-cbind(SM.vals,SM.probs)
    upper.tails<-cbind(rev(SM.dist[,1]),cumsum(rev(SM.dist[,2])))
    outp$cutoff.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),1]
    outp$true.alpha.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),2]  
  }

  
  if(outp$method=="Monte Carlo"){
    mc.perm<-matrix(ncol=outp$k,nrow=outp$n)
    mc.stats<-numeric(n.mc)
    for(i in 1:n.mc){
      for(j in 1:outp$n){
        mc.perm[j,!is.na(outp$obs.mat[j,])]<-sample(possible.ranks[j,!is.na(outp$obs.mat[j,])])
      }
      mc.stats[i]<-round(SM.stat(mc.perm),5)
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