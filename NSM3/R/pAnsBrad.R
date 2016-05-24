pAnsBrad <-function(x,y=NA,g=NA,method=NA,n.mc=10000){
  
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
  outp$m<-length(x)
  outp$n<-length(y)
  outp$ties <- (length(c(x,y)) != length(unique(c(x,y))))
  outp$extra=NULL;
  even<-(outp$m+outp$n+1)%%2
	outp$stat.name<-"Ansari-Bradley C"
	
  
  ##When the user doesn't give us any indication of which method to use, try to pick one.
  if(is.na(method)){
    if(outp$ties){
      if(choose(outp$m+outp$n,outp$n)<=10000){
        method<-"Exact"
      }
      if(choose(outp$m+outp$n,outp$n)>10000){
        method<-"Monte Carlo"
      }
    }
    if(!outp$ties){
      if(outp$m+outp$n<=200){
        method<-"Exact"
      }
      if(outp$m+outp$n>200){
        method<-"Asymptotic"
      }
    }
  }
  #####################################################################
  outp$method<-method
  
	
  if(!outp$ties){
	  if(outp$method=="Monte Carlo"){
		  warning("The exact computation will work for large data without ties, so Exact methods are used rather than Monte Carlo.")
      outp$method="Exact"
	  }
			
	  if(outp$method=="Exact"){
	  	tmp<-ansari.test(y,x,"less",exact=T)
	  	tmp2<-ansari.test(y,x,exact=T)
	  }
	  
    if(outp$method=="Asymptotic"){
	  	tmp<-ansari.test(y,x,"less",exact=F)
 	  	tmp2<-ansari.test(y,x,exact=F)
	  }	

	  outp$obs.stat<-as.numeric(tmp$statistic)
	  outp$p.val<-tmp$p.value
	  outp$two.sided<-tmp2$p.value
  }
  
  if(outp$ties){
    if(outp$method!="Asymptotic"){
      
	    our.data<-rbind(c(x,y),c(rep(1,length(x)),rep(0,length(y))))
	    sorted<-our.data[1,order(our.data[1,]) ]
	    x.labels<-our.data[2,order(our.data[1,]) ]

	    N<-length(sorted)

	    med<-ceiling(N/2)
	    if(N%%2==0){no.ties<-c(1:med,med:1)}
	    if(N%%2==1){no.ties<-c(1:med,(med-1):1)}
	

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

	    assigned.scores<-tied.ranks
 	    outp$obs.stat<-sum(tied.ranks[x.labels==0])
      
    }
    
    if(outp$method=="Exact"){
          possible.orders<-combinations(outp$m+outp$n,outp$n)
          C.stats<-apply(possible.orders,1,function(x) sum(assigned.scores[x]))
          C.tab<-table(C.stats)
          C.vals<-round(as.numeric(names(C.tab)),5)
          C.probs<-as.numeric(C.tab)/sum(C.tab)
          outp$p.val<-sum(C.probs[C.vals>=round(outp$obs.stat,5)])
          outp$two.sided<-2*min(outp$p.val,1-outp$p.val)
    }
    
    
    if(outp$method=="Monte Carlo"){
      outp$n.mc<-n.mc;
      outp$p.val<-0
      for(i in 1:n.mc){
        if(sum(sample(assigned.scores,outp$n))>=outp$obs.stat){
          outp$p.val<-outp$p.val+1/n.mc      
        }
      }
      outp$two.sided<-2*min(outp$p.val,1-outp$p.val)
    }
    
    if(outp$method=="Asymptotic"){
      options(warn = (-1)); tmp<-ansari.test(y,x,"greater",exact=F); options(warn = (0));
      tmp2<-ansari.test(y,x,exact=F)
      outp$obs.stat<-as.numeric(tmp$statistic)
      outp$p.val<-tmp$p.value
      outp$two.sided<-tmp2$p.value
    }	
    
    
  }
  
  class(outp)<-"NSM3Ch5p"
	if(!even&&outp$method=="Exact")
	{outp$extra<-"(N is odd so the null distribution is not symmetric
		and so the two-sided p-value is approximate.)"}
	outp
}
