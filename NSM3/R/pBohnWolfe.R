pBohnWolfe<-function(x,y,k,q,c,d,method="Monte Carlo",n.mc=10000){
  outp<-list()
  outp$m<-length(x)
  outp$n<-length(y)
  outp$n.mc<-n.mc
  if(k*c!=outp$m){warning("Check k*c is the same as the length of x")}
  if(q*d!=outp$n){warning("Check q*d is the same as the length of y")}
  outp$stat.name<-"Bohn-Wolfe U"
  
  outp$method<-method
  
  if(outp$method=="Asymptotic"){warning("The Asymptotic distribution is not yet supported in this version.")}
  
  if(outp$method=="Exact"){warning("The Exact distribution is not yet supported in this version.")}
    outp$method="Monte Carlo"
                           
  mc.dist<-numeric(n.mc)
  
  outp$obs.stat<-0
  for(j in 1:(q*d)){
    outp$obs.stat<-outp$obs.stat+sum(x<y[j])
  }
                             
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
    
  outp$p.val<-sum(mc.probs[mc.vals>=outp$obs.stat])  
      
  class(outp)<-"NSM3Ch5p"  
  outp 
}