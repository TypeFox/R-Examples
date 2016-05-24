pPairedWilcoxon<-function(x,y=NA,g=NA,method=NA,n.mc=10000){
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
  outp$n.mc<-n.mc
  
  if(outp$m!=outp$n){stop("The two groups must be the same size for the paired test.")}
  outp$stat.name<-"Wilcoxon T+"
  
  
  if(!is.na(method)){
    if(method=="Asymptotic"){
      stop("Use wilcox.test() for the asymptotic distribution.")
    }
  }
    
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
  
  z<-y-x
  #Check for zeroes
  num.zero<-sum(z==0)
  if(num.zero>0){
    z<-z[z!=0]
    outp$n<-outp$m<-(outp$m-num.zero)
  }
  psi<-(z>0)
  #Rank will handle tied Z values
  r<-rank(abs(z))
  
  T.calc<-function(psi.tmp){
    return(sum(psi.tmp*r))
  }
  
  
  outp$obs.stat<-T.calc(psi)
  
  if(outp$method=="Exact"){
    possible.comb<-expand.grid(lapply(1:outp$n, function(i) (c(0,1))))
    T.values<-apply(possible.comb,1,T.calc)
    outp$p.val<-mean(T.values>=outp$obs.stat)
  }
  if(outp$method=="Monte Carlo"){
    outp$p.val<-0
    for(i in 1:outp$n.mc){
      outp$p.val<-outp$p.val+1/outp$n.mc*(outp$obs.stat<=T.calc(sample(0:1,outp$n,replace=T)))
    }
  }
  class(outp)<-"NSM3Ch5p"
  outp
}