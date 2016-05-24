sim2.bd.mrca.reconst <-
function(mrca,numbsim,lambda,mu,frac,K){
treearray<- list()
	for (j in 1:numbsim) {
		temp <- 0
		while (temp == 0) {
			phy <- 0
			#while (class(phy[[1]]) != "phylo") {
				phy <- sim2.bd.mrca(mrca,1,lambda,mu,K)
			#}
			phy <- reconstructed.age(phy,frac)
			if (class(phy[[1]])=="phylo"){
				if (max(branching.times.complete(phy[[1]]))>(mrca-0.0001)){   #Tanja 8.7.14 unnecessary!? maybe. or: in the reconstructed tree there may be the mrca deleted
					treearray <- c(treearray,list(phy[[1]]))
					temp = 1
					}
			}		
	}	
	}
	treearray
	}

