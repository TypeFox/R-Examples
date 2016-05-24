"emp.grammar" <-
function(h2object,data,snpsubset,idsubset,strata,times=200,quiet=FALSE,bcast=10) {
	out <- grammar(h2object=h2object,data=data,snpsubset=snpsubset,idsubset=idsubset,strata=strata,
			times=times,quiet=quiet,bcast=bcast)
	out
}

