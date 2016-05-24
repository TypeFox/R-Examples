cLepage <-
function(alpha,m,n, method=NA, n.mc=10000){
	
	outp<-list()
	outp$stat.name<-"Lepage D"
	
  if(alpha>1||alpha<0||class(alpha)!="numeric"){
	  cat('Error: Check alpha value! \n')
	  return(alpha)
	}
  
	outp$alpha<-alpha
	outp$m<-m 
	outp$n<-n 
	outp$n.mc<-n.mc
	N<-m+n
  
	##When the user doesn't give us any indication of which method to use, try to pick one.
	if(is.na(method)){
	  if(choose(N,m)<=10000){
	    method<-"Exact"
	  }
	  if(choose(N,m)>10000){
	    method<-"Monte Carlo"
	  }
	  
	}
	#####################################################################
	outp$method<-method
  
	 
	##Only needs to depend on y values
	D.calc<-function(C.vals,W.vals){
	  odd<-N%%2 
	  if(!odd){
	    exp.C=n*(N+2)/4
	    var.C=m*n*(N+2)*(N-2)/(48*(N-1))
	  }
	  if(odd){
	    exp.C=n*(N+1)^2/(4*N)
	    var.C=m*n*(N+1)*(3+N^2)/(48*N^2)
	  }
	  W.obs<-sum(W.vals)
	  W.star<-(W.obs-n*(N+1)/2)/sqrt(m*n*(N+1)/12)
	  C.star<-(sum(C.vals)-exp.C)/sqrt(var.C)
	  return(W.star^2+C.star^2)
	}
	
  tmp.W<-1:N
	med<-ceiling(N/2)
	if(!N%%2){tmp.C<-c(1:med,med:1)}
	if(N%%2){tmp.C<-c(1:med,(med-1):1)}
	  
	if(outp$method=="Exact"){
	  possible.orders<-combinations(outp$m+outp$n,outp$n)
	  
	  possible.C<-t(apply(possible.orders,1,function(x) tmp.C[x]))
	  possible.W<-t(apply(possible.orders,1,function(x) tmp.W[x]))
    
	  theor.dist<-numeric(nrow(possible.C))
	  for(i in 1:nrow(possible.C)){
	    theor.dist[i]<-D.calc(possible.C[i,],possible.W[i,])
	  }
    
		cutoff.candidates<-sort(unique(theor.dist))
		upper.calc<-function(cand){
			mean(cand<=theor.dist)
		}

		upper.tails<-unlist(lapply(cutoff.candidates,upper.calc))
		outp$cutoff.U<-cutoff.candidates[min(which(upper.tails<=alpha))]
		outp$true.alpha.U<-upper.tails[min(which(upper.tails<=alpha))]
	}

	if(outp$method=="Asymptotic"){
		outp$cutoff.U<-qchisq(1-alpha,2)
	}

	if(outp$method=="Monte Carlo"){
		mc.dist<-numeric(n.mc)
		
    for(i in 1:n.mc){
      mc.sample<-sample(1:N,n)
      mc.dist[i]<-D.calc(tmp.C[mc.sample],tmp.W[mc.sample])
		}
    
		cutoff.candidates<-sort(unique(mc.dist))
		upper.calc<-function(cand){
			mean(cand<=mc.dist)
		}
		
		upper.tails<-unlist(lapply(cutoff.candidates,upper.calc))
		outp$cutoff.U<-cutoff.candidates[min(which(upper.tails<=alpha))]
		outp$true.alpha.U<-upper.tails[min(which(upper.tails<=alpha))]
	}	

	class(outp)="NSM3Ch5c"
	outp
}
