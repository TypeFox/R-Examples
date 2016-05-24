sim2.bd.ssa <-
function(n,numbsim,lambda,mu){
	phy <- list()
	for (j in 1:numbsim){
	stop=0
	while (stop==0) {
		phy2 <- sim2.bd(n,0,lambda,mu)
		if (class(phy2) == "phylo") { 
		phytmp <- drop.extinct(phy2,tol = 0.0000001)
		if (length(phytmp$tip)== n) {
			stop = 1
			}
		}
		}
		phy <- c(phy, list(phy2))
	}
	phy
	}

