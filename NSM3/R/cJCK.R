cJCK <-function(alpha,n, method=NA, n.mc=10000){
    
    outp<-list()
    outp$stat.name<-"Jonckheere-Terpstra J"
    
    
    if(alpha>1||alpha<0||class(alpha)!="numeric"){
      cat('Error: Check alpha value! \n')
      return(alpha)
    }
    
    outp$alpha<-alpha
    outp$n<-n 
    outp$n.mc<-n.mc
    N<-sum(n)
    k<-length(n)
    ##When the user doesn't give us any indication of which method to use, try to pick one.
    if(is.na(method)){
      if(N<=200){
        method<-"Exact"
      }
      if(N>200){
        method<-"Asymptotic"
      }
    }
    #####################################################################   
    
    outp$method<-method
    
    
    if(outp$method=="Monte Carlo"){
      warning("The exact computation will work for large data, so Monte Carlo methods
		    are not recommended for this procedure.")
      outp$method="Exact"
    }
   
    if(outp$method=="Exact"){
      num.comb<-factorial(N)/prod(factorial(outp$n))
      max.J=0;
      for(i in 1:(k-1)){
        max.J=max.J+outp$n[i]*(cumsum(outp$n)[k]-cumsum(outp$n)[i])
      }
      
      #Remember we need to include 0 as possibility, so following code may appear strange at first;
      if(max.J%%2){
        even=1;
        upper=(max.J-1)/2
      }
      
      if(!max.J%%2){
        even=0;
        upper=max.J/2
      }
      
      ##Initialize##
      freq.dist<-numeric(upper+1);
      freq.dist[1]<-1;
      
      
      ##Function##
      update<-function(m,n){
        size.check<-(n+1)<=upper;
        if(size.check){
          p=min(m+n,upper)
          for(t in (n+1):p){
            for(u in upper:t){
              freq.dist[u+1]<<-freq.dist[u+1]-freq.dist[u+1-t]
            }
          }  
        }
        
        q=min(m,upper)
        
        for(s in 1:q){
          for(u in s:upper){
            freq.dist[u+1]<<-freq.dist[u+1]+freq.dist[u+1-s]
          }
        }
      }
      
      for(i in 1:(k-1)){
        update(outp$n[i],(cumsum(outp$n)[k]-cumsum(outp$n)[i]))
      }
      
      low.prob.dist<-freq.dist/num.comb;
      if(even){
        prob.dist<-c(low.prob.dist,rev(low.prob.dist))
      }
      
      if(!even){
        prob.dist<-c(low.prob.dist,rev(low.prob.dist)[-1])
      }
      
      values<-(0:max.J)
      JT.dist<-cbind(values,prob.dist)
          
          
      upper.tails<-cbind(rev(values),cumsum(rev(prob.dist)))
            
      outp$cutoff.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),1]
      
      outp$true.alpha.U<-upper.tails[max(which(upper.tails[,2]<=alpha)),2]
    }
    
    if(outp$method=="Asymptotic"){
      outp$stat.name<-"Jonckheere-Terpstra J*"
      outp$cutoff.U<-qnorm(1-alpha)
    }
        
    class(outp)="NSM3Ch6c"
    outp
  }