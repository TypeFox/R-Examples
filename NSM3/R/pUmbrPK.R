pUmbrPK<-function(x,peak=NA,g=NA,method=NA, n.mc=10000){
  
  if(is.na(peak)){
    warning("The peak is required for this procedure. If peak is unkown, use pUmbrPU instead.")
    return;
  }    
  
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
  outp$stat.name<-paste("Mack-Wolfe Peak Known A",peak)
  
  outp$n.mc<-n.mc
  outp$ties<- (length(x) != length(unique(x)))
  
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
  
  A.calc<-function(obs.data){
    U.vec<-numeric((peak*(peak-1)+(k-peak+1)*(k-peak))/2)
    U.calc<-function(i,j){
      wilcox.test(obs.data[g==levels(g)[i]],obs.data[g==levels(g)[j]])$statistic
    }
  
    count<-0
    for(i in 2:peak){
      for(j in 1:(i-1))	{
        count<-count+1
        options(warn = (-1));
        U.vec[count]<-U.calc(i,j)
        options(warn = (0));
      }
    }
    for(i in peak:(k-1)){
      for(j in (i+1):k)	{
        count<-count+1
        options(warn = (-1));
        U.vec[count]<-U.calc(i,j)
        options(warn = (0));
      }
    }
    sum(U.vec)  
  }
  
  outp$obs.stat<-A.calc(x);
  
  if(!outp$ties){
    if(outp$method=="Monte Carlo"){
      warning("The exact computation will work for large data, so Monte Carlo methods
      		are not recommended for this procedure.")
      outp$method="Exact"
    }
  
    if(outp$method=="Exact"){
      num.comb<-factorial(N)/prod(factorial(outp$n))
      
      max.A=0;
      for(i in 1:(peak-1)){
        max.A=max.A+outp$n[i]*(cumsum(outp$n)[peak]-cumsum(outp$n)[i])
      }
      for(i in peak:(k-1)){
        max.A=max.A+outp$n[i]*(cumsum(outp$n)[k]-cumsum(outp$n)[i])
      }
        
      #Remember we need to include 0 as possibility, so following code may appear strange at first;
      if(max.A%%2){
        even=1;
        upper=(max.A-1)/2
      }
        
      if(!max.A%%2){
        even=0;
        upper=max.A/2
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
        
      update(outp$n[peak],cumsum(outp$n)[k]-outp$n[peak])
        
      if(peak>2){
        for(i in 1:(peak-2)){
          update(outp$n[i],(cumsum(outp$n)[peak-1]-cumsum(outp$n)[i]))
        }
      }
      if(k>peak+1){
        for(i in k:(peak+2)){
          update(outp$n[i],(cumsum(outp$n)[i-1]-cumsum(outp$n)[peak]))
        }
      }
        
      low.freq.dist<-freq.dist;
      if(even){
        freq.dist<-c(low.freq.dist,rev(low.freq.dist))
      }
      if(!even){
        freq.dist<-c(low.freq.dist,rev(low.freq.dist)[-1])
      }
        
      prob.dist<-freq.dist/sum(freq.dist)
      
      values<-(0:max.A)
      A.dist<-cbind(values,prob.dist)
      outp$p.val<-sum(A.dist[values>=outp$obs.stat,2])
    }
    
    
    
    if(outp$method=="Asymptotic"){
      N1<-sum(l[1:peak])
      N2<-sum(l[peak:k])
      A.star<-(outp$obs.stat-(N1^2+N2^2-sum(l^2)-l[peak]^2)/4)/
        sqrt((2*(N1^3+N2^3)+3*(N1^2+N2^2)-sum(l^2*(2*l+3))-l[peak]^2*(2*l[peak]+3)+12*l[peak]*N1*N2-12*l[peak]^2*sum(l))/72)
      outp$stat.name<-paste("Mack-Wolfe Peak Known A*",peak)
      outp$obs.stat<-A.star
      outp$p.val<-1-pnorm(A.star)
    }
  }
  if(outp$ties){
    possible.ranks<-as.numeric(rank(x))
    
    if(outp$method=="Asymptotic"){
      N1<-sum(l[1:peak])
      N2<-sum(l[peak:k])
      A.star<-(outp$obs.stat-(N1^2+N2^2-sum(l^2)-l[peak]^2)/4)/
        sqrt((2*(N1^3+N2^3)+3*(N1^2+N2^2)-sum(l^2*(2*l+3))-l[peak]^2*(2*l[peak]+3)+12*l[peak]*N1*N2-12*l[peak]^2*sum(l))/72)
      outp$stat.name<-paste("Mack-Wolfe Peak Known A*",peak)
      outp$obs.stat<-A.star
      outp$p.val<-1-pnorm(A.star) 
      
      outp$extra<-paste("Ties are present, so exact variance of A",peak," is not available. Reported 
                          asymptotic p-value is therefore larger than the truth.")
    }
    
    if(outp$method=="Exact"){
      if(factorial(N)/prod(factorial(outp$n))>10000){
        use.MC<-(-1)
        while(use.MC<0){
          use.MC <- as.integer(readline("Ties are present, so exact distributions may be computationally intensive. \n Press 0 for Exact or 1 for Monte Carlo computations."))
          use.MC <- ifelse(grepl("\\D",use.MC),NA,as.integer(use.MC))
          if(is.na(use.MC)){break}  # breaks when hit enter
        }
      }
      if(factorial(N)/prod(factorial(outp$n))<=10000){use.MC=0}
      
      if(use.MC){
        outp$method<-"Monte Carlo"
      }
      if(!use.MC){
        possible.orders<-multComb(l)
        possible.perms<-t(apply(possible.orders,1,function(x) possible.ranks[x]))
        
        A.stats<-apply(possible.perms,1,A.calc)
        A.tab<-table(A.stats)
        A.vals<-round(as.numeric(names(A.tab)),5)
        A.probs<-as.numeric(A.tab)/sum(A.tab)
        outp$p.val<-sum(A.probs[A.vals>=round(outp$obs.stat,5)])
      }
    }
    
    if(outp$method=="Monte Carlo"){
      outp$p.val<-0
      for(i in 1:n.mc){
        mc.perm<-sample(possible.ranks)
        if(A.calc(mc.perm)>=outp$obs.stat){
          outp$p.val=outp$p.val+1/n.mc
        }
      }
    }
  }
    
  class(outp)<-"NSM3Ch6p"
  outp
}

