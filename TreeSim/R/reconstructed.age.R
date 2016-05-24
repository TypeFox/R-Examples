reconstructed.age <-
function(treearray,frac){
	phy <- list()
	for (j in 1:length(treearray)){	
		temp <- treearray[[j]]
		if (class(temp)=="phylo"){
		temp<-drop.extinct(temp,tol = 0.00001)
		if (class(temp) != "phylo") {phy <- c(phy, list(1-round(1-frac)))} else {
		ntotal<- length(temp$tip.label)
		delete <- round((1-frac)*ntotal)
		if (ntotal == delete) {
			phy <- c(phy, list(0))
		} else if (ntotal == (delete+1) ) {
			phy <- c(phy, list(1))
		} else {
		temp$tip.label <- paste("t", sample(ntotal), sep = "")
		droptips<-vector()
		if (delete>0){
		for (j in 1:delete){
			droptips<-c(droptips,paste("t", (ntotal+1-j), sep = ""))
			}
		temp<-drop.tip(temp,droptips)
		}
		phy <- c(phy, list(temp))	}
		}} else {
		phy <- c(phy, list(temp))	
		}
	}
	phy
	}

