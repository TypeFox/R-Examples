sim2.bd.fast.age <-
function(n,numbsim,lambda,mu,rho,age,mrca=FALSE){	phy <- list()
	for (j in 1:numbsim){
		if (mrca==FALSE) {
			temp1 <- sim2.bd.fast.single.origin(n,lambda,mu,rho,age)
			temp<-collapse.singles(temp1)
			} else {
			temp <- sim2.bd.fast.single.mrca(n,lambda,mu,rho,age)
			} 
		phy <- c(phy, list(temp))
		}
	phy
}

