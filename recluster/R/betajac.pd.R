betajac.pd <- function(betapd) {
	UniFrac <- (2*betapd$pd.obs.tot-betapd$sum.pd.obs)/betapd$pd.obs.tot
	UniFrac_turn <- 2*(betapd$min.pd.obs)/(betapd$sum.pd.obs-betapd$pd.obs.tot+2*betapd$min.pd.obs)
	UniFrac_PD <- UniFrac-UniFrac_turn
	return(as.matrix(data.frame(UniFrac,UniFrac_turn,UniFrac_PD,row.names=betapd$labcomb)))
}

