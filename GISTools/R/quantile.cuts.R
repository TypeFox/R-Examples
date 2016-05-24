`quantileCuts` <-
function(x,n=5,params=NA) {
	if (!is.na(n)) probs = (1:(n-1))/n	
	if (!is.na(params)) probs = params
	quantile(x,probs)}

