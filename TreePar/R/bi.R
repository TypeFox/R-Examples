bi <- function(i,t,l,mu,psi,rho) {
	if (i==1) {
			out<- (1-2*(1-rho[i]))*l[i]+mu[i]+psi[i]
			out<-out/ai(i,l,mu,psi)	
	} else {
			out<- (1-2*(1-rho[i])*p(i-1,t[i],t,l,mu,psi,rho))*l[i]+mu[i]+psi[i]
			out<-out/ai(i,l,mu,psi)	}
	#print("B")
	#print(c(out,i))
	out
	}
