get.weights.unstab <- function(x, stop.method = NULL, withSampW = TRUE){
	if(!((class(x) == "iptw") | class(x) == "mniptw")) stop("\"get.weights.unstab\" is only defined for iptw objects.")
	if(class(x) == "iptw"){
	prodWt <- rep(1, length(x$psList[[1]]$treat))
	for(i in 1:length(x$psList)){
		hdWt <- x$psList[[i]]$ps * x$psList[[i]]$treat + (1-x$psList[[i]]$ps) * (1-x$psList[[i]]$treat)
		prodWt <- hdWt * prodWt
	}	
	return(1/prodWt)
	}
	
	if(class(x) == "mniptw"){
		denom <- 1/get.weights(x$psList[[1]], stop.method = stop.method, estimand = "ATE", withSampW = FALSE)
		for(i in 2:length(x$psList)) denom <- denom * 1/get.weights(x$psList[[i]], stop.method = stop.method, estimand = "ATE", withSampW = FALSE)
		w <- 1/denom
		if(withSampW) w <- w * x$psList[[1]]$sampw
	}
}