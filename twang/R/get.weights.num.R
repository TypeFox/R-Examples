get.weights.num <- function(iptw, fitList){
	hdWt <- rep(1, length(iptw$psList[[1]]$treat))
	for(i in 1:length(iptw$psList)){
		hdPred <- fitted(fitList[[i]]) * iptw$psList[[i]]$treat + (1-fitted(fitList[[i]])) * (1-iptw$psList[[i]]$treat)
		hdWt <- hdWt * hdPred
	}
	return(hdWt)
}