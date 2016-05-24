pHaySton<-function(x,g=NA,method=NA,n.mc=10000){
  outp<-list()
  outp$stat.name<-"Hayter-Stone W*"
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
  outp$n<-l
  
  outp$num.comp<-num.comp<-k*(k-1)/2	
  outp$ties<- (length(x) != length(unique(x))) 
  
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
  for(i in 1:(k-1)){
    for(j in (i+1):k){
      outp$labels[count]<-paste(levels(g)[i],"-",levels(g)[j])
      count<-count+1
    }
  }
  
  W.star.calc<-function(x,i,j){
    group.sizes<-l[c(i,j)]
    W.stat<-sum(rank(c(x[g==levels(g)[i]],x[g==levels(g)[j]]))[(group.sizes[1]+1):sum(group.sizes)])  
    W.mean<-group.sizes[2]*(sum(group.sizes)+1)/2
    tie.vec<-as.numeric(table(c(x[g==levels(g)[i]],x[g==levels(g)[j]])))
    W.var<-prod(group.sizes)/24*(sum(group.sizes)+1-sum((tie.vec-1)*tie.vec*(tie.vec+1)/(sum(group.sizes)*(sum(group.sizes)-1))))
    (W.stat-W.mean)/sqrt(W.var)
  }
  
  W.star.all<-function(x){
    W.star.vec<-numeric(num.comp)
    count<-1
    for(i in 1:(k-1)){
      for(j in (i+1):k){
        W.star.vec[count]<-W.star.calc(x,i,j)
        count<-count+1
      }
    }
    W.star.vec
  }
  
  outp$obs.stat<-W.star.all(x);
  
  if(outp$method=="Exact"){
    possible.combs<-multComb(l)
    if(outp$ties){
      possible.ranks<-as.numeric(rank(x))
      possible.combs<-t(apply(possible.combs,1,function(x) possible.ranks[x]))
    }
    exact.dist<-apply(possible.combs,1,function(x) max(W.star.all(x)))
    for(i in 1:num.comp){
      outp$p.val[i]<-mean(exact.dist>=outp$obs.stat[i])
    }
  }
  if(outp$method=="Asymptotic"){
	if(length(unique(outp$n))==1){    
		for(j in 1:num.comp){
      		outp$p.val[j]<-pHayStonLSA((outp$obs.stat[j]),k)
    		}
	}
  	if(length(unique(outp$n))!=1){
		warning("Since sample sizes are unequal, Hayter-Stone LSA is not appropriate.")
		outp$method="Monte Carlo"
	}
  }  
  
  if(outp$method=="Monte Carlo"){
    outp$p.val<-numeric(num.comp)
    for(i in 1:n.mc){
      mc.order<-as.numeric(sample(rank(x)))
      for(j in 1:num.comp){
        if(max(W.star.all(mc.order))>=outp$obs.stat[j]){
          outp$p.val[j]=outp$p.val[j]+1/n.mc
        }
      }
    }
  }  

  class(outp)<-"NSM3Ch6MCp"
  outp  
}

