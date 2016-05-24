cNDWol<-function(alpha,n, method=NA, n.mc=10000){
  outp<-list()
  outp$stat.name<-"Nemenyi, Damico-Wolfe Y*"
  
  if(alpha>1||alpha<0||class(alpha)!="numeric"){
    cat('Error: Check alpha value! \n')
    return(alpha)
  }    
  
  
  outp$n.mc<-n.mc  
  
  l<-n
  k<-length(n)
  N<-sum(n)
  outp$trt<-n[1]
  outp$n<-n[-1]
  g<-rep(1:k,n)
  
  outp$num.comp<-num.comp<-k-1
    
  ##When the user doesn't give us any indication of which method to use, try to pick one.
  if(is.na(method)){
    if(factorial(sum(n))/prod(factorial(n))<=10000){
      method<-"Exact"
    }
    if(factorial(sum(n))/prod(factorial(n))>10000){
      method<-"Monte Carlo"
    }
  }
  #####################################################################
  
  outp$method<-method
      
  gcd<- function(u, v) {
    a<-max(u,v)
    b<-min(u,v)
    for(i in 1:b){
      if(a%%(b/i)==0){
        return(b/i)
      }
    }
  }
  
  #Note that someone has already used the name "lcm"
  LCM<-function(u,v){
    u*v/gcd(u,v)
  }
  N.star<-LCM(outp$trt,outp$n[1])
  if(k>2){
    for(i in 2:(k-1)){
      N.star<-LCM(N.star,outp$n[i])
    }
  }
  
  possible.ranks<-1:N
  
  Y.star.calc<-function(obs.order){
    R.vec<-unlist(lapply(1:k,function(x) mean(obs.order[g==x])))
    N.star*(R.vec[-1]-R.vec[1])
  }
  
  if(outp$method=="Exact"){
    possible.combs<-multComb(l)
    Y.stats<-round(apply(possible.combs,1,function(x) max(Y.star.calc(x))),10)
    Y.tab<-table(Y.stats)
    Y.vals<-as.numeric(names(Y.tab))
    Y.probs<-as.numeric(Y.tab)/sum(Y.tab)
    Y.dist<-cbind(Y.vals,Y.probs)
    #From cJCK.R;
    upper.tails<-cbind(rev(Y.vals),cumsum(rev(Y.probs)))
    outp$cutoff.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),1]
    outp$true.alpha.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),2]
  }
  
  if(outp$method=="Monte Carlo"){
    mc.dist<-numeric(outp$n.mc)
    for(i in 1:outp$n.mc){
      mc.dist[i]<-round(max(Y.star.calc(sample(1:N))),5)
    }
    mc.values<-sort(unique(mc.dist))
    mc.probs<-as.numeric(table(mc.dist))/outp$n.mc
    Y.dist<-cbind(mc.values,mc.probs)
    
    #From cJCK.R;
    upper.tails<-cbind(rev(mc.values),cumsum(rev(mc.probs)))
    outp$cutoff.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),1]
    outp$true.alpha.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),2]    
  }
  if(outp$method=="Asymptotic"){
    if(length(unique(outp$n))==1){
      outp$cutoff.U<-cMaxCorrNor(alpha,k-1,outp$n[1]/(outp$trt+outp$n[1]))*sqrt(N*(N+1)/12)*sqrt(1/outp$trt+1/outp$n[1])
    }
    if(length(unique(outp$n))>1){
      outp$cutoff.U<-qnorm(1-alpha/(k-1))*sqrt(N*(N+1)/12)*sqrt(1/outp$trt+1/min(outp$n))
    }
    
  }    
  class(outp)<-"NSM3Ch6MCc"
  outp
}


