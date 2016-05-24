RFPW<-function(z){
	n<-length(z)

	f.star<-function(tmp,i,j,k){
		sign(tmp[i]+tmp[j]-2*tmp[k])+sign(tmp[i]+tmp[k]-2*tmp[j])+sign(tmp[j]+tmp[k]-2*tmp[i])
	}
	
	#Compute the T statistic;
	f.vals<-numeric(n*(n-1)*(n-2)/6)
	iter<-1
	for(i in 1:(n-2)){
		for(j in (i+1):(n-1)){
			for(k in (j+1):n){
				f.vals[iter]<-f.star(z,i,j,k)
				iter<-iter+1						
			}
		}
	}
	T.stat<-sum(f.vals)

	#Compute the univariate B statistics;
	b.calc<-function(s){
		A<-B<-E<-0;
		if(s<(n-1)){
			for(j in (s+1):(n-1)){
				for(k in (j+1):n){
					A=A+f.star(z,s,j,k)
				}
			}
		}
		if(s>1&&s<n){		
			for(j in 1:(s-1)){
				for(k in (s+1):n){
					B<-B+f.star(z,j,s,k)
				}
			}
		}
		if(s>2){	
			for(j in 1:(s-2)){
				for(k in (j+1):(s-1)){
					E<-E+f.star(z,j,k,s)
				}
			}
		}
	return(A+B+E)
	}	
	b.vals<-unlist(lapply(1:n,b.calc))
	
	#Compute the bivariate B statistics;
	b2.calc<-function(u,v){
		A<-B<-E<-0;
		if(u>1){
			for(j in 1:(u-1)){
					A<-A+f.star(z,j,u,v)
			}
		}
		if((v-u)>1){		
			for(j in (u+1):(v-1)){
					B<-B+f.star(z,u,j,v)
			}
		}
		if(v<n){	
			for(j in (v+1):n){
					E<-E+f.star(z,u,v,j)
			}
		}
	return(A+B+E)
	}
	b2.vals<-numeric(n*(n-1)/2)
	iter<-1
	for(i in 1:(n-1)){
		for(j in (i+1):n){
			b2.vals[iter]<-b2.calc(i,j)
			iter<-iter+1
		}
	}

	sigma2.hat<-(n-3)*(n-4)/((n-1)*(n-2))*sum(b.vals^2)+(n-3)/(n-4)*sum(b2.vals^2)+n*(n-1)*(n-2)/6-(1-(n-3)*(n-4)*(n-5)/(n*(n-1)*(n-2)))*T.stat^2

  outp<-list()
  outp$obs.stat<-T.stat/sqrt(sigma2.hat)
  outp$p.val<-2*(1-pnorm(abs(outp$obs.stat)))
  outp
}
