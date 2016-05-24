Fderifuncshift<-function(time,t,lambda,mu,rho=1) {  #page 117 in TPB paper
	k <- inter(time,t)
	out1 <- Fderifuncshifth(time,t,lambda,mu,k)
	out2 <- (1-rho+rho*(Ffuncshifth(time,t,lambda,mu)+1))^2
	out<-rho*out1/out2
	out
}