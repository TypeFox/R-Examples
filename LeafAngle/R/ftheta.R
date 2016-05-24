# Vectorized version of fthetaf
ftheta <- function(distpars=NA, ...){
	res <- mapply(fthetaf, ..., MoreArgs=list(distpars=distpars))
return(res)
}