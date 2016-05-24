pFligPoli <-function(x,y=NA,g=NA,method=NA,n.mc=10000){
  
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
  outp$stat.name<-"Fligner-Policello U"
  
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

  possible.ranks<-rank(c(x,y))

  U.calc<-function(rank.vec){
    x.tmp<-rank.vec[1:outp$m]
    y.tmp<-rank.vec[(outp$m+1):(outp$m+outp$n)]
    p.vec<-unlist(lapply(x.tmp,function(x){sum(x>y.tmp)+0.5*sum(x==y.tmp)}))
    q.vec<-unlist(lapply(y.tmp,function(y){sum(y>x.tmp)+0.5*sum(y==x.tmp)}))
    p.bar<-mean(p.vec)
    q.bar<-mean(q.vec)
    v1<-sum((p.vec-p.bar)^2)
    v2<-sum((q.vec-q.bar)^2)
    return((sum(q.vec)-sum(p.vec))/(2*sqrt(v1+v2+p.bar*q.bar)))
  }
  
  outp$obs.stat<-U.calc(possible.ranks)
  
  if(outp$method=="Exact"){
          possible.orders<-multComb(c(outp$m,outp$n))
          possible.perm<-t(apply(possible.orders,1,function(x){possible.ranks[x]}))
          U.stats<-apply(possible.perm,1,U.calc)
          U.tab<-table(U.stats)
          U.vals<-round(as.numeric(names(U.tab)),5)
          U.probs<-as.numeric(U.tab)/sum(U.tab)
          outp$p.val<-sum(U.probs[U.vals>=round(outp$obs.stat,5)])
          outp$two.sided<-2*min(outp$p.val,1-outp$p.val)
    }
    
    if(outp$method=="Monte Carlo"){
      outp$n.mc<-n.mc;
      outp$p.val<-0
      for(i in 1:n.mc){
        if(U.calc(sample(possible.ranks))>=outp$obs.stat){
          outp$p.val<-outp$p.val+1/n.mc      
        }
      }
      outp$two.sided<-2*min(outp$p.val,1-outp$p.val)
    }
    
    if(outp$method=="Asymptotic"){
      outp$p.val<-1-pnorm(outp$obs.stat)
      outp$two.sided<-2*min(outp$p.val,1-outp$p.val)
    }	
  
  class(outp)<-"NSM3Ch5p"
	outp
}
