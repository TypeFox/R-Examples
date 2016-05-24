getExpectedValues <-
function(x){
#http://sites.stat.psu.edu/online/courses/stat504/R/death.R
	if( !is.array(x) ){ warning(paste0(sQuote("x"), " must be matrix or array")) }
	lendm <- length(dim(x))

	margins <- getMargins(x)
	if( any(unlist(margins) == 0) ){ margins <- lapply(margins, function(v){v+1}); }
	sm <- sum(margins[[1]])

	E <- switch(as.character(lendm)
		, '2' = (margins[[1]]%o%margins[[2]])/(sm)
		, '3' = (margins[[1]]%o%margins[[2]]%o%margins[[3]])/(sm^2)
		, eval(parse(text=paste("margins[[", 1:lendm, "]]", collapse="%o%")))/(sm^(lendm-1))
		)
	dimnames(E) <- dimnames(x)
return(E)
}
