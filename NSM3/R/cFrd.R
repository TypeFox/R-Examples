cFrd<-function(alpha, k, n, method=NA, n.mc=10000){
  outp<-list()
  outp$stat.name<-"Friedman, Kendall-Babington Smith S"
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
  
  S.calc<-function(x){
    sum.squares<-sum(colSums(t(apply(x,1,rank)))^2)
    return(round(12/(n*k*(k+1))*sum.squares-3*n*(k+1),10))  
  }
  
  possible.ranks<-matrix(rep(1:outp$k,outp$n),ncol=outp$k,byrow=T)
  
  if(outp$method=="Exact"){
    
      phi<-function(full){
        mat<-full[,-1]
        sort<-t(apply(mat,1,sort))
        uniq<-unique(sort,MARGIN=1)
        sort<-cbind(full[,1],sort)
        counts<-numeric(nrow(uniq))
        for(i in 1:length(counts)){
          counts[i]<-sum(apply(sort,1,function(x,y) x[1]*identical(x[-1],y), y=uniq[i,]))
        }
        return(as.matrix(cbind(counts[],uniq)))
      }
      
      update<-function(full,original.with.ranks){
        mat<-full[,-1]
        original<-original.with.ranks[,-1]
        output<-matrix(nrow=dim(original)[1]*max(nrow(full),1),ncol=min(dim(mat)[2],length(mat))+1)
        count<-1
        for(i in 1:max(dim(mat)[1],1)){
          if(max(dim(mat)[1],1)==1){
            for(j in 1:max(dim(original)[1],1)){
              output[count,]<-c(full[1],mat+original[j,])
              count<-count+1
            }
          }
          if(max(dim(mat)[1],1)!=1){
            for(j in 1:max(dim(original)[1],1)){
              output[count,]<-c(full[i,1],mat[i,]+original[j,])  
              count<-count+1
            }
          }
        }
        return(output)
      }
      
      
      initial<-cbind(rep(1,factorial(k)),multComb(rep(1,k)))
      if(nrow(initial)!=factorial(k)){
        return("Error!")
      }
      
      new<-update(phi(initial),initial)
      if(n>2){
        for(i in 1:(n-2)){
          new<-phi(update(new,initial))
        }
      }
      sum.squares<-apply(new[,-1]^2,1,sum)
      S.vals<-round(12/(n*k*(k+1))*sum.squares-3*n*(k+1),5)
      S.probs<-new[,1]/sum(new[,1])
      
      S.vals2<-unique(S.vals)
      S.probs2<-numeric(length(S.vals2))
      for(i in 1:length(S.vals2)){
        S.probs2[i]<-sum(S.probs[S.vals==S.vals2[i]])
      }
      S.dist<-cbind(sort(S.vals2),S.probs2[order(S.vals2)])
      upper.tails<-cbind(rev(S.dist[,1]),cumsum(rev(S.dist[,2])))
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
    mc.stats[i]<-S.calc(mc.perm)
  }
  
  mc.vals<-sort(unique(mc.stats))
  mc.dist<-as.numeric(table(mc.stats))/n.mc
  
  upper.tails<-cbind(rev(mc.vals),cumsum(rev(mc.dist)))
  outp$cutoff.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),1]
  outp$true.alpha.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),2]
  
}

if(outp$method=="Asymptotic"){
  outp$cutoff.U<-qchisq(1-alpha,outp$k-1)
}

class(outp)<-"NSM3Ch7c"
  return(outp)
}