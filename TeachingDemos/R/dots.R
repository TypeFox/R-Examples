"dots" <-
function(x,...){
	
	sx <- sort(x)
	sy <- unlist(lapply(table(sx),seq))
	plot(sx,sy, xlab=deparse(substitute(x)), ylab="Count",...)
	
}

