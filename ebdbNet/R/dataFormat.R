`dataFormat` <- 
function(longitudinal.data) 
{
	if(length(unique(attributes(longitudinal.data)$repeats)) > 1) {
		stop("Error: number of replicates must remain constant over time.")
	}
	P <- attributes(longitudinal.data)$dim[2]
	T <- length(attributes(longitudinal.data)$time)
	R <- attributes(longitudinal.data)$repeats[1]

	dat <- vector("list", R)
	for(r in 1:R) {
		dat[[r]] <- t(longitudinal.data[seq(from = r, to = R*T, by = R),])
	}

	return(dat)
}