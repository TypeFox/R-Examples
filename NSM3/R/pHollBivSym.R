pHollBivSym<-function(x,y=NA,g=NA,method=NA,n.mc=10000){
	
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
  
  if(outp$m!=outp$n){stop("'x' and 'y' must have the same number of observations")}
  
  outp$n.mc<-n.mc 
  outp$stat.name<-"Hollander A"

  if(!is.na(method)){
  if(method=="Asymptotic"){
    warning("Koziol's LSA was found to perform poorly so Asymptotic method not included.")
    outp$method=NA
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
  
  
	obs.data<-cbind(x,y)
	a.vec<-apply(obs.data,1,min)
	b.vec<-apply(obs.data,1,max)

	test<-function(r,c) {as.numeric((a.vec[c]<b.vec[r])&&(b.vec[r]<=b.vec[c])&&(a.vec[r]<=a.vec[c]))}
	myVecFun <- Vectorize(test,vectorize.args = c('r','c')) 
	d.mat<-outer(1:outp$n, 1:outp$n, FUN=myVecFun) 

	A.calc<-function(r.vec){
		s.vec<-2*r.vec-1
		T.vec<-s.vec%*%d.mat
		A.obs<-sum(T.vec*T.vec)/n^2
		return(A.obs)
	}

	outp$obs.stat<-A.calc(apply(obs.data,1,function(x){x[1]<x[2]}))

if(outp$method=="Exact"){
	possible.r<-expand.grid(lapply(1:n, function(i) (c(0,1)))) 
	A.values<-apply(possible.r,1,A.calc)
	outp$p.val<-mean(A.values>=outp$obs.stat)
}

if(outp$method=="Monte Carlo"){
  outp$p.val<-0
	for(i in 1:n.mc){
		mc.r<-sample(0:1,n,replace=T)
		outp$p.val<-outp$p.val+1/n.mc*(A.calc(mc.r)>=outp$obs.stat)
	}	
}  

  class(outp)="NSM3Ch5p"
  outp

}

