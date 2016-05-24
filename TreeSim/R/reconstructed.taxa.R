reconstructed.taxa <-
function(treearray,delete){
	phy <- list()
	for (j in 1:length(treearray)){
		temp <- reconstructed.taxa.single(treearray[[j]],delete)
		phy <- c(phy, list(temp))
		}
	phy
	}

