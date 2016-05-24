pNDWol<-function(x,g=NA,method=NA, n.mc=10000){
  outp<-list()
  outp$stat.name<-"Nemenyi, Damico-Wolfe Y*"
  outp$n.mc<-n.mc  
  
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
  
  outp$trt<-l[1]
  outp$n<-l[-1]
  
  
  outp$num.comp<-num.comp<-k-1
  outp$ties <- (length(x) != length(unique(x)))
    
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

  count<-1
  outp$labels<-character(num.comp)
  for(j in 2:k){
    outp$labels[count]<-paste("1-",levels(g)[j])
    count<-count+1
  }
  
  
  gcd<- function(u, v) {
    a<-max(u,v)
    b<-min(u,v)
    for(i in 1:b){
      if(a%%(b/i)==0){
        return(b/i)
      }
    }
  }
  
  #Note that someone has already used the name "lcm"
  LCM<-function(u,v){
    u*v/gcd(u,v)
  }
  N.star<-LCM(outp$trt,outp$n[1])
  if(k>2){
    for(i in 2:(k-1)){
      N.star<-LCM(N.star,outp$n[i])
    }
  }
  possible.ranks<-rank(x)
  
  Y.star.calc<-function(obs.order){
    R.vec<-unlist(lapply(1:k,function(x) mean(obs.order[g==levels(g)[x]])))
    N.star*(R.vec[-1]-R.vec[1])
  }
  outp$obs.stat<-Y.star.calc(possible.ranks)
    
  if(outp$method=="Exact"){
    possible.combs<-multComb(l)
    if(outp$ties){
      possible.combs<-t(apply(possible.combs,1,function(x) possible.ranks[x]))
    }
    exact.dist<-apply(possible.combs,1,function(x) max(Y.star.calc(x)))
    for(i in 1:num.comp){
      outp$p.val[i]<-mean(exact.dist>=outp$obs.stat[i])
    }
  }
  
  if(outp$method=="Monte Carlo"){
    mc.dist<-numeric(outp$n.mc)
    for(i in 1:outp$n.mc){
      mc.dist[i]<-max(Y.star.calc(sample(1:N)))
    }
    for(i in 1:num.comp){
      outp$p.val[i]<-mean(mc.dist>=outp$obs.stat[i])
    }
  }
  if(outp$method=="Asymptotic"){
    if(length(unique(outp$n))==1){
      for(i in 1:num.comp){
        outp$p.val[i]<-pMaxCorrNor(outp$obs.stat[i]/(N.star*sqrt(N*(N+1)/12)*sqrt(1/outp$trt+1/outp$n[i])),k-1,outp$n[1]/(outp$trt+outp$n[1]))
      }
    }
    if(length(unique(outp$n))>1){
      for(i in 1:num.comp){
        outp$p.val[i]<-pnorm(outp$obs.stat[i]/(N.star*sqrt(N*(N+1)/12)*sqrt(1/outp$trt+1/outp$n[i])),k-1,outp$n[1]/(outp$trt+outp$n[1]),lower.tail=F)/(k-1)
      }
    }
    
  }    
  class(outp)<-"NSM3Ch6MCp"
  outp
}


