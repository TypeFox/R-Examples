checkPars<-function(x,n) {
if ((length(x) > 1) && (length(x) != n))
	stop(paste("length of",deparse(substitute(x)),"vector does not match number of variables !"))
if (length(x) == 1) return(rep(x,n))
return(x)
}