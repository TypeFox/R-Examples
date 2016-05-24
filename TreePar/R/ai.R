ai <-function(i,l,mu,psi){
	temp<-(l[i]-mu[i]-psi[i])^2
	out<-sqrt(temp + 4 * l[i]*psi[i])
	#print("A")
	#print(out)
	out
	}
