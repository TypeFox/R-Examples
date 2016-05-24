reconstructed.taxa.single <-
function(tree,delete){
		temp <- tree
		temp<-drop.extinct(temp,tol = 0.000001)
		ntotal<- length(temp$tip.label)
		if (ntotal == delete) {
			phy = 0
		} else if (ntotal == (delete+1) ) {
			phy = 1
		} else {
			temp$tip.label <- paste("t", sample(length(temp$tip)), sep = "")
			droptips<-vector()
			if (delete>0){
			for (j in 1:delete){
				droptips<-c(droptips,paste("t", (ntotal+1-j), sep = ""))
			}
			temp<-drop.tip(temp,droptips)
			}
			phy <- temp
		}
	phy
	}

