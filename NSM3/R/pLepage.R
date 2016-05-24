pLepage <-
  function(x,y=NA,g=NA,method=NA,n.mc=10000){
    
    ##Adapted from kruskal.test()##
    if (is.list(x)) {
      if (length(x) < 2L) 
        stop("'x' must be a list with at least 2 elements")
      y<-x[[2]]
      x<-x[[1]]
    }
    else {
      if(min(is.na(y))!=0){
        k<-length(unique(g))
      if (length(x) != length(g)) 
        stop("'x' and 'g' must have the same length")
      if (k < 2) 
        stop("all observations are in the same group")
      y<-x[g==2]
      x<-x[g==1]
      }
    }
    #####################
        
    
    outp<-list()
    outp$m<-m<-length(x)
    outp$n<-n<-length(y)
    outp$n.mc<-n.mc
    N<-outp$m+outp$n
    outp$ties <- (length(c(x,y)) != length(unique(c(x,y))))
    even<-(outp$m+outp$n+1)%%2
    outp$stat.name<-"Lepage D"

    
    ##When the user doesn't give us any indication of which method to use, try to pick one.
    if(is.na(method)){
      if(choose(outp$m+outp$n,outp$n)<=10000){
        method<-"Exact"
      }
      if(choose(outp$m+outp$n,outp$n)>10000){
        method<-"Monte Carlo"
      }
    }
    #####################################################################
    outp$method<-method
        
    
    tmp.W<-rank(c(x,y))
  
    our.data<-rbind(c(x,y),c(rep(1,length(x)),rep(0,length(y))))
  	sorted<-our.data[1,order(our.data[1,]) ]
  	x.labels<-our.data[2,order(our.data[1,]) ]
	  	
  	med<-ceiling(N/2)
  	if(even){no.ties<-c(1:med,med:1)}
  	if(!even){no.ties<-c(1:med,(med-1):1)}
	
	
  	obs.group<-numeric(N)
  	group.num<-1
  	for(i in 1:N){
  	  if(obs.group[i]==0){
  	    obs.group[i]<-group.num
  	    for(j in i:N){
  	      if(sorted[i]==sorted[j]){
  	        obs.group[j]<-obs.group[i]
  	      }
  	    }
  	    group.num<-group.num+1;
  	  }
  	}
	
  	group.ranks<-tapply(no.ties,obs.group,mean)
	
  	tied.ranks<-numeric(N)
  	for(i in 1:group.num){
  	  tied.ranks[which(obs.group==as.numeric(names(group.ranks)[i]))]<-group.ranks[i]  
  	}
	
  tmp.C<-c(tied.ranks[x.labels==1],tied.ranks[x.labels==0])
    
  
  ##Only needs to depend on y values
	D.calc<-function(C.vals,W.vals){
		 
		if(even){
			exp.C=n*(N+2)/4
			var.C=m*n*(N+2)*(N-2)/(48*(N-1))
		}
		if(!even){
			exp.C=n*(N+1)^2/(4*N)
			var.C=m*n*(N+1)*(3+N^2)/(48*N^2)
		}
		W.obs<-sum(W.vals)
		W.star<-(W.obs-n*(N+1)/2)/sqrt(m*n*(N+1)/12)
		C.star<-(sum(C.vals)-exp.C)/sqrt(var.C)
		return(W.star^2+C.star^2)
	}

	outp$obs.stat<-D.calc(tmp.C[(m+1):N],tmp.W[(m+1):N])

	if(outp$method=="Exact"){
	  possible.orders<-combinations(outp$m+outp$n,outp$n)
	  
    possible.C<-t(apply(possible.orders,1,function(x) tmp.C[x]))
	  possible.W<-t(apply(possible.orders,1,function(x) tmp.W[x]))
    
    theor.dist<-numeric(nrow(possible.C))
    for(i in 1:nrow(possible.C)){
		  theor.dist[i]<-D.calc(possible.C[i,],possible.W[i,])
    }
    
		outp$p.val<-mean(theor.dist>=outp$obs.stat)
  }

	if(outp$method=="Asymptotic"){
		outp$p.val<-(1-pchisq(outp$obs.stat,2))
	}
	
  
  if(outp$method=="Monte Carlo"){
		outp$p.val<-0
		for(i in 1:n.mc){
      mc.sample<-sample(1:N,n)
      
			if(D.calc(tmp.C[mc.sample],tmp.W[mc.sample])>=outp$obs.stat){
				outp$p.val=outp$p.val+1/n.mc
			}
		}
	}	
	class(outp)="NSM3Ch5p"
	outp
}
