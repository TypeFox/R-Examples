"emp.qtscore" <-
function(formula,data,snpsubset,idsubset,strata,trait.type="gaussian",times=200,quiet=FALSE,bcast=10) {
	out <- qtscore(formula=formula,data=data,snpsubset=snpsubset,idsubset=idsubset,strata=strata,
			trait.type=trait.type,times=times,quiet=quiet,bcast=bcast)
	out
}

