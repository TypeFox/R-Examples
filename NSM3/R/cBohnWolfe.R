cBohnWolfe<-function(alpha,k,q,c,d,method="Monte Carlo",n.mc=10000){
  if(alpha>1||alpha<0||class(alpha)!="numeric"){
    cat('Error: Check alpha value! \n')
    return(alpha)
  }
  outp<-list()
  outp$m<-k*c
  outp$n<-q*d
  outp$n.mc<-n.mc
  outp$stat.name<-"Bohn-Wolfe U"
  outp$method<-method
  
  if(outp$method=="Asymptotic"){warning("The Asymptotic distribution is not yet supported in this version.")}
  
  if(outp$method=="Exact"){warning("The Exact distribution is not yet supported in this version.")}
  
  outp$method="Monte Carlo"
                           
  mc.dist<-numeric(n.mc)
                           
  for(iter in 1:n.mc){
    sample<-NULL
    for(j in 1:c){
      for(i in 1:(k)){
        sample<-c(sample,rbeta(1,i,k+1-i))
      }
    }
    for(j in 1:d){
      for(i in 1:(q)){
        sample<-c(sample,rbeta(1,i,q+1-i))
      }  
    } 
    stat<-0
    for(j in (k*c+1):(k*c+q*d)){
      stat<-stat+sum(sample[1:(k*c)]<sample[j])
    }
    mc.dist[iter]<-stat
  }
  mc.vals<-as.numeric(names(table(mc.dist)))
  mc.probs<-table(mc.dist)/n.mc
  upper.tails<-cbind(rev(mc.vals),cumsum(rev(mc.probs)))
  
  outp$cutoff.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),1]
  outp$true.alpha.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),2]  
  
  class(outp)<-"NSM3Ch5c"
  outp 
}