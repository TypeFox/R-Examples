pJCK<-function(x,g=NA,method=NA, n.mc=10000){
   
  ##From kruskal.test()##
  if (is.list(x)) {
    if (length(x) < 2L) 
      stop("'x' must be a list with at least 2 elements")
    DNAME <- deparse(substitute(x))
    x <- lapply(x, function(u) u <- u[complete.cases(u)])
    k <- length(x)
    l <- sapply(x, "length")
    if (any(l == 0)) 
      stop("all groups must contain data")
    g <- factor(rep(1:k, l))
    x <- unlist(x)
  }
  else {
    if (length(x) != length(g)) 
      stop("'x' and 'g' must have the same length")
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
    OK <- complete.cases(x, g)
    x <- x[OK]
    g <- g[OK]
    if (!all(is.finite(g))) 
      stop("all group levels must be finite")
    g <- factor(g)
    l<-as.numeric(table(g))
    k <- nlevels(g)
    if (k < 2) 
      stop("all observations are in the same group")
  }
  N <- length(x)
  #####################
  outp<-list()
  outp$n<-l
  outp$stat.name<-"Jonckheere-Terpstra J"
  
  outp$n.mc<-n.mc
  outp$ties <- (length(x) != length(unique(x)))
  
  ##When the user doesn't give us any indication of which method to use, try to pick one.
  if(is.na(method)){
    if(outp$ties){
      if(factorial(sum(outp$n))/prod(factorial(outp$n))<=10000){
        method<-"Exact"
      }
      if(factorial(sum(outp$n))/prod(factorial(outp$n))>10000){
        method<-"Monte Carlo"
      }
    }
    if(!outp$ties){
      method<-"Exact"
    }
  }
  #####################################################################
  
  outp$method<-method
  
  JT.calc<-function(obs.data){
    U.vec<-numeric(k*(k-1)/2)
    U.calc<-function(i,j){
      wilcox.test(obs.data[g==levels(g)[i]],obs.data[g==levels(g)[j]])$statistic
    }
    
    count<-0
    for(i in 2:k){
      for(j in 1:(i-1))	{
        count<-count+1
        options(warn = (-1));
        U.vec[count]<-U.calc(i,j)
        options(warn = (0));
      }
    }
    sum(U.vec)
  }
  
  outp$obs.stat<-JT.calc(x);
  
  if(!outp$ties){
    
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
      outp$p.val<-sum(JT.dist[values>=outp$obs.stat,2])
    }
    
    if(outp$method=="Asymptotic"){
      J.star<-(outp$obs.stat-(N^2-sum(l^2))/4)/sqrt((N^2*(2*N+3)-sum(l^2*(2*l+3)))/72)
      outp$stat.name<-"Jonckheere-Terpstra J*"
      outp$obs.stat<-J.star
      outp$p.val<-1-pnorm(J.star)
    }
  }
  
  if(outp$ties){  
    possible.ranks<-as.numeric(rank(x))
    
    if(outp$method=="Asymptotic"){
      tie.counts<-as.numeric(table(possible.ranks))
      J.var<-1/72*(N*(N-1)*(2*N+5)-sum(outp$n*(outp$n-1)*(2*outp$n+5))-
                     sum(tie.counts*(tie.counts-1)*(2*tie.counts+5)))+
             1/(36*N*(N-1)*(N-2))*sum(outp$n*(outp$n-1)*(outp$n-2))*
                     sum(tie.counts*(tie.counts-1)*(tie.counts-2))+
             1/(8*N*(N-1))*sum(outp$n*(outp$n-1))*sum(tie.counts*(tie.counts-1))
              
      J.star<-(outp$obs.stat-(N^2-sum(l^2))/4)/sqrt(J.var)
      outp$stat.name<-"Jonckheere-Terpstra J*"
      outp$obs.stat<-J.star
      outp$p.val<-1-pnorm(J.star)
    }
      
    if(outp$method=="Exact"){
        possible.orders<-multComb(l)
        possible.perms<-t(apply(possible.orders,1,function(x) possible.ranks[x]))
                  
        J.stats<-apply(possible.perms,1,JT.calc)
        J.tab<-table(J.stats)
        J.vals<-round(as.numeric(names(J.tab)),5)
        J.probs<-as.numeric(J.tab)/sum(J.tab)
        outp$p.val<-sum(J.probs[J.vals>=round(outp$obs.stat,5)])
        
      }
    if(outp$method=="Monte Carlo"){
      outp$p.val<-0
      for(i in 1:n.mc){
        mc.perm<-sample(possible.ranks)
        if(JT.calc(mc.perm)>=outp$obs.stat){
          outp$p.val=outp$p.val+1/n.mc
        }
      }
    }
  }
  
  class(outp)<-"NSM3Ch6p"
  outp
}

