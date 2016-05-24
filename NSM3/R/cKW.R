cKW <-
  function(alpha,n, method=NA, n.mc=10000){
    
    outp<-list()
    outp$stat.name<-"Kruskal-Wallis H"
    
    if(alpha>1||alpha<0||class(alpha)!="numeric"){
      cat('Error: Check alpha value! \n')
      return(alpha)
    }
        
    outp$alpha<-alpha
    outp$n<-n 
    outp$n.mc<-n.mc
    N<-sum(n)
    
    ##When the user doesn't give us any indication of which method to use, try to pick one.
    if(is.na(method)){
      if(factorial(N)/prod(factorial(outp$n))<=10000){
        method<-"Exact"
      }
      if(factorial(N)/prod(factorial(outp$n))>10000){
        method<-"Monte Carlo"
      }
    }
    #####################################################################   
        
    outp$method<-method
    
    H.calc<-function(obs.data,k){
      N<-sum(k)
      tmp<-cumsum(k)
      tmp2<-cumsum(rank(obs.data))[tmp]
      R.2<-c(tmp2[1]^2,diff(tmp2)^2)/k
      12/(N^2+N)*sum(R.2)-3*N-3  
    }
    
    
    if(outp$method=="Exact"){
      possible.comb<-multComb(n)
      theor.dist<-round(apply(possible.comb,1,H.calc,k=n),8)
      cutoff.candidates<-sort(unique(theor.dist))
      upper.calc<-function(cand){
        mean(cand<=theor.dist)
      }
      
      upper.tails<-unlist(lapply(cutoff.candidates,upper.calc))
      outp$cutoff.U<-cutoff.candidates[min(which(upper.tails<=alpha))]
      outp$true.alpha.U<-upper.tails[min(which(upper.tails<=alpha))]
    }
    
    if(outp$method=="Asymptotic"){
      outp$cutoff.U<-qchisq(1-alpha,length(n)-1)
    }
    
    if(outp$method=="Monte Carlo"){
      mc.dist<-numeric(n.mc)
      for(i in 1:n.mc){
        mc.dist[i]<-round(H.calc(sample(1:N),n),8)
      }
      cutoff.candidates<-sort(unique(mc.dist))
      
      upper.calc<-function(cand){
        mean(cand<=mc.dist)
      }
      
      upper.tails<-unlist(lapply(cutoff.candidates,upper.calc))
      outp$cutoff.U<-cutoff.candidates[min(which(upper.tails<=alpha))]
      outp$true.alpha.U<-upper.tails[min(which(upper.tails<=alpha))]
      
    }	
    
    
    
    class(outp)="NSM3Ch6c"
    outp
  }
