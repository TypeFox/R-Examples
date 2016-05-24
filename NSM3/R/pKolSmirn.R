pKolSmirn <-
function(x,y=NA,g=NA,method=NA,n.mc=10000){
    
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
  outp$n.mc<-n.mc
  N<-outp$m+outp$n
  outp$stat.name<-"Kolmogorov-Smirnov J"
  
  
  ##When the user doesn't give us any indication of which method to use, try to pick one.
  if(is.na(method)){
    if(outp$m+outp$n<=1000){
      method<-"Exact"
    }
    if(outp$m+outp$n>1000){
      method<-"Asymptotic"
    }
    
  }
  #####################################################################
  outp$method<-method 

	gcd<- function(u, v) {
		a<-max(u,v)
		b<-min(u,v)
		for(i in 1:b){
			if(a%%(b/i)==0){
				return(b/i)
			}
		}
	}

	d<-gcd(outp$m,outp$n)
	
  tmp.ks1<-ks.test(x,y)
	outp$obs.stat<-as.numeric(tmp.ks1$statistic)*d

	if(outp$method=="Monte Carlo"){
		warning("The exact computation will work for large data, so Monte Carlo methods
				are not recommended for this procedure.")
	  outp$method<-"Exact"
  }

	if(outp$method=="Exact"){
		outp$p.val<-tmp.ks1$p.val
	}

	if(outp$method=="Asymptotic"){
		outp$stat.name<-"Kolmogorov-Smirnov J*"
		outp$obs.stat<-d*outp$obs.stat/(outp$m*outp$n*N)^(1/2)
		Qfun<-function(s){
			k<-c(-100:100)
			q<-(-1)^k*exp(-2*k^2*s^2)
			round(1-sum(q),4)
		}
		outp$p.val<-(Qfun(outp$obs.stat))
	}
	class(outp)<-"NSM3Ch5p"
	outp
}
