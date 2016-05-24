makeDeltaRanksTable <-
function(ll) {
	dats=names(ll)
	kognames=ll[[1]][,1]
	kogtable=c()
	for (n in dats){
		kogtable=data.frame(cbind(kogtable,ll[[n]][,3][match(kognames,ll[[n]][,1])]))
	}
	names(kogtable)=dats
	row.names(kogtable)=kognames
	return(kogtable)
}
