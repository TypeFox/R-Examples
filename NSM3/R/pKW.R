pKW<-function(x,g=NA, method=NA, n.mc=10000){
  
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
  outp$stat.name<-"Kruskal-Wallis H"
  outp$n.mc<-n.mc
  outp$obs.stat<-kruskal.test(x,g)$statistic
    
  
  ##When the user doesn't give us any indication of which method to use, try to pick one.
  if(is.na(method)){
      if(factorial(sum(outp$n))/prod(factorial(outp$n))<=10000){
        method<-"Exact"
      }
      if(factorial(sum(outp$n))/prod(factorial(outp$n))>10000){
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
  
  outp$ties <- (length(x) != length(unique(x)))
  
  if(outp$ties){
    possible.ranks<-sort(as.numeric(rank(x)))
    if(method=="Asymptotic"){
        outp$stat.name<-"Kruskal-Wallis H'"
    }
    if(method!="Asymptotic"){
      outp$obs.stat=H.calc(x,l)
    }
  }  
  
  if(method=="Asymptotic"){
    outp$p.val<-kruskal.test(x,g)$p.val  
  }
  
  

    
  if(method=="Exact"){
    size<-factorial(sum(l))/prod(factorial(l))
    if(!outp$ties){
      possible.perms<-multComb(l)
    }
    if(outp$ties){
      possible.orders<-multComb(l)
      possible.perms<-t(apply(possible.orders,1,function(x) possible.ranks[x]))
    }
    H.stats<-apply(possible.perms,1,H.calc,k=l)
    outp$p.val<-mean(H.stats>=outp$obs.stat)
  }

  if(method=="Monte Carlo"){
    outp$p.val<-0
    for(i in 1:n.mc){
      if(!outp$ties){
        mc.perm<-sample(1:N)
      }
      if(outp$ties){
        mc.perm<-sample(possible.ranks)
      }
      if(H.calc(mc.perm,l)>=outp$obs.stat){
        outp$p.val=outp$p.val+1/n.mc
      }
    }
  }
  class(outp)<-"NSM3Ch6p"
  outp
}

