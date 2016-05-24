mra.control <- function(algorithm=1, maxfn=1000, cov.meth=1, trace=0, tol=1e-7 ){  #, rel.tol=1e-10 ){

if( !is.numeric(algorithm) || !(algorithm %in% c(1)) ){
	warning("value of 'algorithm' should be 1 (VA09AD).  Other values ignored.")
	algorithm <- 1
}

if( !is.numeric(maxfn) || maxfn <=0 ){
 	stop("maximum number of function evaluations must be > 0")
}

if( !is.numeric(cov.meth) || !(cov.meth %in% c(1,2)) ){
 	stop("value of 'cov.meth' must be either 1 or 2")
}

if( !is.numeric(trace) ){
	stop("trace must be an numeric")
}

if( !is.numeric(tol) ){
    stop("tolerance must be numeric")
}

list( algorithm=algorithm[1], maxfn=round(maxfn[1]), 
	cov.meth=cov.meth[1], trace=trunc(trace[1]), tol=tol)  
}
