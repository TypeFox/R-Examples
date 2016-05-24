sim2.bd.rateshift <- function(n,numbsim,lambda,mu,frac,times,K,norm){	phy <- list()
	time<-vector()
	for (j in 1:numbsim){
		temp <- sim2.bd.rateshift.single(n,lambda,mu,frac,times,K,norm)
		phy <- c(phy, list(temp[[1]]))
		time<-c(time,temp[[2]])
		}
	phy2<-list(phy,time)
	phy2
}
