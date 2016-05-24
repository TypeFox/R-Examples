sim2.bd.fast.single <- function(n,lambda,mu,rho){
	r <- runif(1,0,1)
	if (lambda>mu) {
	#origin <- log((-lambda * rho - lambda * r^(1/n) + mu * r^(1/n) + lambda * rho * r^(1/n))/(lambda * rho * (-1 + r^(1/n)))) / (lambda - mu)
	#bug 30.4.14
	lambda2<-lambda*rho
	mu2<-mu-lambda*(1-rho)
	origin <- log((lambda2-mu2 * r^(1/n))/(lambda2*(1- r^(1/n))))/(lambda2-mu2)
	} else {
	#origin<- -(r^(1/n)/(lambda *(-1 + r^(1/n)* rho))) %bug 16.4.14
	origin<- -(r^(1/n)/(lambda * rho*(-1 + r^(1/n))))
	}
	phy<-sim2.bd.fast.single.origin(n,lambda,mu,rho,origin)
	phy2<-collapse.singles(phy)
	phy2<-list(phy2,origin)
	phy2
	}
