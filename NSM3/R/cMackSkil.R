cMackSkil<-function(alpha,k,n,c, method=NA, n.mc=10000){
  outp<-list()
  outp$stat.name<-"Mack-Skillings MS"
  outp$n.mc<-n.mc
  
  outp$k<-k
  outp$n<-n
  outp$c<-c
  
  if(alpha>1||alpha<0||class(alpha)!="numeric"){
    cat('Error: Check alpha value! \n')
    return(alpha)
  }
  
  outp$alpha<-alpha
  
  num.obs<-outp$k*outp$n*outp$c
  
  ##When the user doesn't give us any indication of which method to use, try to pick one.
  if(is.na(method)){
    if(factorial(outp$c*outp$k)^outp$n<=10000){
      method<-"Exact"
    }
    if(factorial(outp$c*outp$k)^outp$n>10000){
      method<-"Monte Carlo"
    }
  }
  #####################################################################
  
  outp$method<-method
  
  possible.ranks<-t(matrix(rep(1:(outp$c*outp$k),outp$n),ncol=outp$n,byrow=F))
  
  MS.calc<-function(obs.data){
    S.vec<-NULL
    for(i in 1:outp$k){
      S.vec<-c(S.vec,sum(obs.data[,((i-1)*outp$c+1):(i*outp$c)])/outp$c)
    }
    MS.stat<-12/(outp$k*(num.obs+outp$n))*sum(S.vec^2)-3*(num.obs+outp$n)
    return(MS.stat)
  }
  
  
  if(outp$method=="Exact"){
    possible.perm<-multCh7(possible.ranks)
    exact.dist<-apply(possible.perm,3,MS.calc)
    
    MS.vals<-sort(unique(exact.dist))
    MS.probs<-as.numeric(table(exact.dist))/(factorial(outp$c*outp$k)^outp$n)
    MS.dist<-cbind(MS.vals,MS.probs)
    upper.tails<-cbind(rev(MS.dist[,1]),cumsum(rev(MS.dist[,2])))
    outp$cutoff.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),1]
    outp$true.alpha.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),2]  
  }
  
  if(outp$method=="Monte Carlo"){
    mc.perm<-matrix(ncol=(outp$c*outp$k),nrow=outp$n)
    mc.stats<-numeric(n.mc)
    for(i in 1:n.mc){
      for(j in 1:outp$n){
        mc.perm[j,]<-sample(possible.ranks[j,])
      }
      mc.stats[i]<-round(MS.calc(mc.perm),5)
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
  