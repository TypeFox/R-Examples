sim.rateshift.taxa.loop <- function(n,numbsim,lambda,mu,frac,times,complete=TRUE,K,norm){
	phy <- sim2.bd.rateshift(n,numbsim,lambda,mu,frac,times,K,norm)
	phy2<-phy[[1]]
	if (complete == FALSE) {
		phy2 <- reconstructed.taxa(phy[[1]],(round(n/frac[1])-n))
		}
	phy2<-list(phy2,phy[[2]])
	phy2
	}
