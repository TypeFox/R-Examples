"mps.ave" <- function(x, method="soerensen", all=FALSE, foc=NULL, what="mean", ...){
	WHAT <- c("mean", "sd")
	what <- pmatch(what, WHAT)
	mat.sim <- sim(x, method=method, ...)
	vek.sim <- as.numeric(mat.sim)
	mean.sim <- mean(vek.sim)
	sd.sim <- sd(vek.sim)
	if(!is.null(foc)){
		if(is.character(foc)){
			foc <- which(rownames(x)==foc)
		}
		mean.sim <- mean(as.matrix(mat.sim)[foc,-foc])
		sd.sim <- sd(as.matrix(mat.sim)[foc,-foc])
	}
	res <- c(mean.sim, sd.sim)
	names(res) <- c("mean", "sd")
	if(!all){
		res <- res[what]
	}
	return(res)
}