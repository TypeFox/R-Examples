sim2.bd.fast <- function(n,numbsim,lambda,mu,rho){	phy <- list()
	time<-vector()
	for (j in 1:numbsim){
		temp <- sim2.bd.fast.single(n,lambda,mu,rho)
		phy <- c(phy, list(temp[[1]]))
		time<-c(time,temp[[2]])
		}
	phy<-list(phy,time)
	phy
}
