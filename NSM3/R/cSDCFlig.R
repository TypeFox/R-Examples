cSDCFlig<-function(alpha,n,method=NA,n.mc=10000){  
  outp<-list()
  outp$stat.name<-"Dwass, Steel, Critchlow-Fligner W"
  
  
  if(alpha>1||alpha<0||class(alpha)!="numeric"){
    cat('Error: Check alpha value! \n')
    return(alpha)
  } 
  
  
  
  outp$n.mc<-n.mc
  outp$n<-l<-n
  outp$alpha<-alpha
  k<-length(n)
  N<-sum(n)
  g<-rep(1:k,n)
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
  
  outp$num.comp<-num.comp<-k*(k-1)/2	
  
  W.star.calc<-function(x,i,j){
    group.sizes<-l[c(i,j)]
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
  
  
  
  if(method=="Exact"){
    possible.combs<-multComb(outp$n)
    exact.dist<-round(apply(possible.combs,1,function(x) max(abs(W.star.all(x)))),9)
    values<-sort(unique(exact.dist))
    prob.dist<-as.numeric(table(exact.dist))/sum(as.numeric(table(exact.dist)))
    
    upper.tails<-cbind(rev(values),cumsum(rev(prob.dist)))
    outp$cutoff.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),1]
    outp$true.alpha.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),2]
  }
  
  
  if(method=="Monte Carlo"){
    mc.dist<-numeric(n.mc)
    for(i in 1:n.mc){
      mc.order<-as.numeric(sample(1:N,N))
      mc.dist[i]<-max(abs(W.star.all(mc.order)))      
    }
    values<-sort(unique(mc.dist))
    prob.dist<-as.numeric(table(mc.dist))/n.mc
    
    upper.tails<-cbind(rev(values),cumsum(rev(prob.dist)))
    outp$cutoff.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),1]
    outp$true.alpha.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),2]
  }  
  if(method=="Asymptotic"){
    outp$cutoff.U<-cRangeNor(alpha,k)
  }
  
  class(outp)<-"NSM3Ch6MCc"
  outp  
}
