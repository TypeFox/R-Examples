sim.gsa.taxa <-
function(treearray,n,frac=1,sampling=1,complete=TRUE){
	if (complete==TRUE) {
		phy<-sim2.gsa(treearray,n,sampling)
		} else  {
		phy<-sim2.gsa(treearray, round(n/frac), sampling)
		for (j in 1:length(phy)){
			phy <- reconstructed.taxa(phy,(round(n/frac)-n))
		}
		}
	phy
	}

