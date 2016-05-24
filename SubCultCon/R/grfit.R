grfit <-
function(dvec,matches,m,l,n){
	grd=1:n
	cl=1-1/l
	for(i in 1:n){
		svec=(m-matches[i,])*dvec/(1-dvec[i]*dvec)-matches[i,]*dvec*cl/(dvec[i]*dvec*cl+1/l)
		grd[i]=sum(svec)-svec[i]
	}
	grd
}
