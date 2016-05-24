plotstemsections <- function(plant, diammultiplier=1, color="chocolate4"){

	stems <- plant$stems
	branches <- plant$branches
	stemlen <- sapply(stems, function(x)sqrt(sum(x$xyz$from - x$xyz$to)^2))
	branchlen <- sapply(branches, function(x)sqrt(sum(x$xyz$from - x$xyz$to)^2))

	# Which sections have non-zero length?
	istemnonzero <- which(stemlen > 0)
	ibranchnonzero <- which(branchlen > 0)

	for(i in istemnonzero){
		st <- stems[[i]]
		if(st$diam > 0){
			plot3dcylinder(start=st$xyz$from, end=st$xyz$to, radius=diammultiplier*st$diam/2)
		}
	}
	for(i in ibranchnonzero){
		st <- branches[[i]]
		if(st$diam > 0){
			plot3dcylinder(start=st$xyz$from, end=st$xyz$to, radius=diammultiplier*st$diam/2)
		}
	}

}
