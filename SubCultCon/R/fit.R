fit <-
function(dvec,matches,m,l,n){
	pimat=1/l+(1-1/l)*(dvec%*%t(dvec))
	smat=matches*log(pimat)+(m-matches)*log(1-pimat)
	val=(sum(diag(smat))-sum(smat))/2
	val
}
