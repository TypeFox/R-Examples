"emp.ccfast" <-
function(y,data,snpsubset,idsubset,times=200,quiet=FALSE,bcast=10) {
	out <- ccfast(y=y,data=data,snpsubset=snpsubset,idsubset=idsubset,
			times=times,quiet=quiet,bcast=bcast)
	out
}

