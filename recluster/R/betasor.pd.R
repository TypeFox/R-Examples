betasor.pd <- function(betapd) {
	PhyloSor <- (2*betapd$pd.obs.tot-betapd$sum.pd.obs)/betapd$sum.pd.obs
	PhyloSor_turn <- betapd$min.pd.obs/(betapd$sum.pd.obs-betapd$pd.obs.tot+betapd$min.pd.obs)
	PhyloSor_PD <- (abs(betapd$dif.pd.obs)/betapd$sum.pd.obs)*((betapd$sum.pd.obs-betapd$pd.obs.tot)/(betapd$sum.pd.obs-betapd$pd.obs.tot+betapd$min.pd.obs))
	return(as.matrix(data.frame(PhyloSor,PhyloSor_turn,PhyloSor_PD,row.names=betapd$labcomb)))
}
