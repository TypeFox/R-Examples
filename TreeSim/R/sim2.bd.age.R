sim2.bd.age <-
function(age,numbsim,lambda,mu,K){
	phy <- list()
	for (j in 1:numbsim){
		temp <- sim2.bd(0,age,lambda,mu,K)
		phy <- c(phy, list(temp))
		}
	phy
	}

