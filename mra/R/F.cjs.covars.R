"F.cjs.covars" <-
function(nan, ns){
	x <- rep( diag(ns), rep(nan, ns*ns) )
	dim(x) <- c(nan, ns, ns)
	dimnames(x) <- list( NULL, paste("Trap",1:ns), paste("X",1:ns, sep="") )
	
	ans <- list( x=x )
	ans
}

