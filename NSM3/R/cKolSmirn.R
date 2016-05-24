cKolSmirn <-
function(alpha,m,n,method=NA,n.mc=10000){
	N=m+n;

	outp<-list()
	outp$m<-m
	outp$n<-n
	outp$alpha<-alpha
	outp$stat.name<-"Kolmogorov-Smirnov J"
	
	##When the user doesn't give us any indication of which method to use, try to pick one.
	if(is.na(method)){
	  if(outp$m+outp$n<=200){
	    method<-"Exact"
	  }
	  if(outp$m+outp$n>200){
	    method<-"Asymptotic"
	  }
	  
	}
	#####################################################################
  outp$method<-method
	
	if(alpha>1||alpha<0||class(alpha)!="numeric"){
	  cat('Error: Check alpha value! \n')
	  return(alpha)
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

	d<-gcd(m,n)

	if(outp$method=="Exact"){
		w.upper<-c(1:outp$m, (outp$m+1):(outp$m+outp$n))
		z.upper<-cumsum(ifelse(order(w.upper) <= outp$m, 1/outp$m, -1/outp$n))
		upper.J<-outp$m*outp$n/d*max(abs(z.upper))
		lower.J<-floor(abs(outp$m-outp$n)/2)
		upper.tail<-rep(0,(upper.J-lower.J+1))
		index<-c(lower.J:upper.J)
		for(i in 1:length(upper.tail)){
			upper.tail[i]<-1 - .Call("pSmirnov2x", p = as.double(index[i]*d/(outp$m*outp$n)), 
			                         as.integer(outp$m), as.integer(outp$n),PACKAGE="NSM3")
		}
		outp$true.alpha.U<-upper.tail[min(which(upper.tail<=alpha))]
		outp$cutoff.U<-index[max(which(upper.tail==outp$true.alpha.U))]
	}

	if(outp$method=="Monte Carlo"){
		warning("The exact computation will work for large data, so Monte Carlo methods
			are not recommended for this procedure.")
	}

	if(outp$method=="Asymptotic"){
		Qfun<-function(s){
			k<-c(-100:100)
			q<-(-1)^k*exp(-2*k^2*s^2)
			round(1-sum(q),4)
		}
		qalpha<-0
		upper<-0
		test<-seq(0.3, 2.39, 0.01)
		for(i in 1:length(test)){
			qalpha[i]<-test[i]
			upper[i]<-Qfun(test[i])
		}
		qtable<-matrix(c(qalpha,upper),ncol=2)
		pos<-which.min(abs(qtable[,2]-alpha))
		J_cut<-qtable[pos,1]
		outp$cutoff.U<-J_cut*(m*n*N)^(1/2)/d
	}

	class(outp)<-"NSM3Ch5c"
	outp
}
