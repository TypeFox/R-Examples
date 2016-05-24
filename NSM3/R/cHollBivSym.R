cHollBivSym<-function(alpha,d.mat,method=NA, n.mc=10000){
	
	outp<-list()
		
  if(alpha>1||alpha<0||class(alpha)!="numeric"){
	  cat('Error: Check alpha value! \n')
	  return(alpha)
	}
  
	outp$alpha<-alpha
	outp$m<-outp$n<-m<-n<-nrow(d.mat)
  	  
  outp$n.mc<-n.mc
	outp$stat.name<-"Hollander A"
  N<-outp$m+outp$n
    
	if(!is.na(method)){
	if(method=="Asymptotic"){
	  warning("Koziol's LSA was found to perform poorly so Asymptotic method not included.")
	  outp$method=NA
	}
	}
	##When the user doesn't give us any indication of which method to use, try to pick one.
	if(is.na(method)){
	  if(choose(N,outp$m)<=10000){
	    method<-"Exact"
	  }
	  if(choose(N,outp$m)>10000){
	    method<-"Monte Carlo"
	  }
	  
	}
	#####################################################################
	outp$method<-method
  
	 
	A.calc<-function(r.vec){
	  s.vec<-2*r.vec-1
	  T.vec<-s.vec%*%d.mat
	  A.obs<-sum(T.vec*T.vec)/n^2
	  return(A.obs)
	}
	  
	if(outp$method=="Exact"){
	  possible.r<-expand.grid(lapply(1:n, function(i) (c(0,1)))) 
	  A.values<-apply(possible.r,1,A.calc)
      
    A.dist<-sort(unique(A.values))
		
    upper.calc<-function(cand){
			mean(cand<=A.values)
		}

		upper.tails<-unlist(lapply(A.dist,upper.calc))
		outp$cutoff.U<-A.dist[min(which(upper.tails<=alpha))]
		outp$true.alpha.U<-upper.tails[min(which(upper.tails<=alpha))]
	}

	if(outp$method=="Monte Carlo"){
		mc.dist<-numeric(n.mc)
		
    for(i in 1:n.mc){
      mc.r<-sample(0:1,outp$n,replace=T)
      mc.dist[i]<-A.calc(mc.r)
		}
    
		A.dist<-sort(unique(mc.dist))
		upper.calc<-function(cand){
			mean(cand<=mc.dist)
		}
		
		upper.tails<-unlist(lapply(A.dist,upper.calc))
		outp$cutoff.U<-A.dist[min(which(upper.tails<=alpha))]
		outp$true.alpha.U<-upper.tails[min(which(upper.tails<=alpha))]
	}	

	class(outp)="NSM3Ch5c"
	outp
}
