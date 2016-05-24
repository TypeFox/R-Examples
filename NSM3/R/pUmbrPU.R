pUmbrPU<-function(x,g=NA,method=NA, n.mc=10000){
  
  
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
  outp$stat.name<-paste("Mack-Wolfe Peak Unknown A*(p-hat)")
  
  outp$n.mc<-n.mc
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
  
  cumulative.sizes<-cumsum(outp$n)
  
  peak.picker<-function(obs.data){
    tmp<-numeric(k)
    for(i in 1:k){
      first<-obs.data[max(1,cumulative.sizes[i-1]+1):cumulative.sizes[i]]
      second<-obs.data[-(max(1,cumulative.sizes[i-1]+1):cumulative.sizes[i])]
      options(warn = (-1));  
      tmp[i]<-(wilcox.test(first,second)$statistic-(outp$n[i]*(cumulative.sizes[k]-outp$n[i])/2))/sqrt(
        outp$n[i]*(cumulative.sizes[k]-outp$n[i])*(cumulative.sizes[k]+1)/12)
      options(warn = (0));
    }
    initial.peak<-peak<-which.max(tmp)
    
    #check for multiple peaks;
    if(1-sum(tmp[-peak]==tmp[peak])){
      return(peak)
    }
    for(i in (initial.peak+1):k){
      if(tmp[i]==tmp[initial.peak]){
        peak<-c(peak,i)
      }
    }
    return(peak)
  }
  
  A.star.calc<-function(obs.data,peak){
    N1<-cumulative.sizes[peak]
    N2<-(cumulative.sizes[k]-max(0,cumulative.sizes[peak-1]))
    exp.Ap<-(N1^2+N2^2-sum(outp$n^2)-outp$n[peak]^2)/4
    var.Ap<-1/72*(2*(N1^3+N2^3)+3*(N1^2+N2^2)-sum(outp$n^2*(2*outp$n+3))
                  -outp$n[peak]^2*(2*outp$n[peak]+3)+12*outp$n[peak]*N1*N2-12*outp$n[peak]^2*cumulative.sizes[k])
    
    
    U.vec<-numeric((peak*(peak-1)+(k-peak+1)*(k-peak))/2)
    U.calc<-function(i,j){
      wilcox.test(obs.data[g==levels(g)[i]],obs.data[g==levels(g)[j]])$statistic
    }
    
    count<-0
    if(peak>1){
    for(i in 2:peak){
      for(j in 1:(i-1))  {
        count<-count+1
        options(warn = (-1));
        U.vec[count]<-U.calc(i,j)
        options(warn = (0));
      }
    }
    }
    if(peak<k){
    for(i in peak:(k-1)){
      for(j in (i+1):k)	{
        count<-count+1
        options(warn = (-1));
        U.vec[count]<-U.calc(i,j)
        options(warn = (0));
      }
    }
    }
    (sum(U.vec)-exp.Ap)/sqrt(var.Ap)
    
  }
  
  PU.calc<-function(obs.data){
    tmp.peak<-peak.picker(obs.data)
    num.peak<-length(tmp.peak)
    tmp.stat<-numeric(num.peak)
    if(num.peak==1){
      tmp.stat<-A.star.calc(obs.data,tmp.peak)
    }
    if(num.peak>1){
      for(i in 1:num.peak){
        tmp.stat[i]<-A.star.calc(obs.data,tmp.peak[i])
      }
    }
    mean(tmp.stat)
  }
  
  
  
  
  est.peak<-peak.picker(x)
  outp$extra<-paste("Estimated Peak Group(s):",est.peak)
  outp$obs.stat<-PU.calc(x)

  possible.ranks<-as.numeric(rank(x))
  
  if(outp$method=="Asymptotic"){
    warning("The asymptotic distribution for this statistic is unknown!")
    outp$method=="Monte Carlo"
  }
  if(outp$method=="Exact"){
    possible.orders<-multComb(l)
    possible.perms<-t(apply(possible.orders,1,function(x) possible.ranks[x]))
    
    A.stats<-apply(possible.perms,1,PU.calc)
    A.tab<-table(A.stats)
    A.vals<-round(as.numeric(names(A.tab)),5)
    A.probs<-as.numeric(A.tab)/sum(A.tab)
    outp$p.val<-sum(A.probs[A.vals>=round(outp$obs.stat,5)])    
  }

  if(outp$method=="Monte Carlo"){
    outp$p.val<-0
    for(i in 1:n.mc){
      mc.perm<-sample(possible.ranks)
      if(PU.calc(mc.perm)>=outp$obs.stat){
        outp$p.val=outp$p.val+1/n.mc
      }
    }
  }

  class(outp)<-"NSM3Ch6p"
  outp
}