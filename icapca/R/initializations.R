initializations <- function(subgaussians, supergaussians, gaussians){
	subs=as.integer(subgaussians)
	if(is.na(subs) || !is.finite(subs) || subs<0){
		stop("unable to coerce subgaussians to a finite non-negative integer")
	}
	supers=as.integer(supergaussians)
	if(is.na(supers) || !is.finite(supers) || supers<0){
		stop("unable to coerce supergaussians to a finite non-negative integer")
	}
	gauss=as.integer(gaussians)
	if(is.na(gauss) || !is.finite(gauss) || gauss<0){
		stop("unable to coerce gaussians to a finite non-ngegative integer")
	}
	m=subs+supers+gauss
	if(m==0) return(0)
	nongaussians=subs+supers
	count=choose(m, nongaussians)*choose(nongaussians, subs)
	return(count)
}