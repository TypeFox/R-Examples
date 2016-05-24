genperms.custom <-
function(numiter=10000,randfun=randfun.default,...) {
	arglist <- list(...)
	return(replicate(numiter,do.call(randfun,arglist)))
	}
