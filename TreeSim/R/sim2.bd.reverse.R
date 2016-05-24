sim2.bd.reverse <- function(n,numbsim,lambda,mu){	phy <- list()
	time<-vector()
	frac<-1
	for (j in 1:numbsim){
		temp <- sim2.bd.reverse.single(n,lambda,mu,frac)
		phy <- c(phy, list(temp[[1]]))
		time<-c(time,temp[[2]])
		}
	phy2<-list(phy,time)
	phy2
}
