cHaySton<-function(alpha,n, method=NA, n.mc=10000){
  outp<-list()
  outp$stat.name<-"Hayter-Stone W*"
  
  
  if(alpha>1||alpha<0||class(alpha)!="numeric"){
    cat('Error: Check alpha value! \n')
    return(alpha)
  }     
  
  outp$n.mc<-n.mc
  
  k<-length(n)
  N<-sum(n)
  outp$n<-n
  g<-rep(1:k,outp$n)   
  
  outp$num.comp<-num.comp<-k*(k-1)/2  
  
  ##When the user doesn't give us any indication of which method to use, try to pick one.
  if(is.na(method)){
    if(factorial(sum(outp$n))/prod(factorial(outp$n))<=10000){
      method<-"Exact"
    }
    if(factorial(sum(outp$n))/prod(factorial(outp$n))>10000){
      method<-"Monte Carlo"
    }
  }
  #####################################################################
  outp$method<-method 
  
  W.star.calc<-function(x,i,j){
    group.sizes<-n[c(i,j)]
    W.stat<-sum(rank(c(x[g==i],x[g==j]))[(group.sizes[1]+1):sum(group.sizes)])  
    W.mean<-group.sizes[2]*(sum(group.sizes)+1)/2
    tie.vec<-as.numeric(table(c(x[g==i],x[g==j])))
    W.var<-prod(group.sizes)/24*(sum(group.sizes)+1-sum((tie.vec-1)*tie.vec*(tie.vec+1)/(sum(group.sizes)*(sum(group.sizes)-1))))
    (W.stat-W.mean)/sqrt(W.var)
  }
  
  W.star.all<-function(x){
    W.star.vec<-numeric(num.comp)
    count<-1
    for(i in 1:(k-1)){
      for(j in (i+1):k){
        W.star.vec[count]<-W.star.calc(x,i,j)
        count<-count+1
      }
    }
    W.star.vec
  }
  
  possible.ranks<-1:N
  
  if(outp$method=="Exact"){
    possible.combs<-multComb(outp$n)
    C.stats<-round(apply(possible.combs,1,function(x) max(W.star.all(x))),10)
    C.tab<-table(C.stats)
    C.vals<-round(as.numeric(names(C.tab)),5)
    C.probs<-as.numeric(C.tab)/sum(C.tab)
    C.dist<-cbind(C.vals,C.probs)
    #From cJCK.R;
    upper.tails<-cbind(rev(C.vals),cumsum(rev(C.probs)))
    outp$cutoff.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),1]
    outp$true.alpha.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),2]   
  }
  if(outp$method=="Asymptotic"){
    if(length(unique(outp$n))==1){    
      test.grid<-seq(0,5,.01)
      for(i in 1:length(test.grid)){
        tmp<-1-pHayStonLSA((test.grid[i]),k,delta=.01)
        if(tmp<=alpha){
          outp$cut.off.U<-test.grid[i]
          outp$true.alpha.U<-tmp
          break
        }
      }
    }
    if(length(unique(outp$n))!=1){
      warning("Since sample sizes are unequal, Hayter-Stone LSA is not appropriate.")
      outp$method="Monte Carlo"
    }
  }  
  
  if(outp$method=="Monte Carlo"){
    mc.dist<-numeric(outp$n.mc)
    for(i in 1:outp$n.mc){
      mc.perm<-sample(possible.ranks)
      mc.dist[i]<-round(max(W.star.all(mc.perm)),10)
    }
    mc.values<-sort(unique(mc.dist))
    mc.probs<-as.numeric(table(mc.dist))/outp$n.mc
    C.dist<-cbind(mc.values,mc.probs)
    
    #From cJCK.R;
    upper.tails<-cbind(rev(mc.values),cumsum(rev(mc.probs)))
    outp$cutoff.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),1]
    outp$true.alpha.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),2]
  
  }  
  
  
  class(outp)<-"NSM3Ch6MCc"
  outp  
}