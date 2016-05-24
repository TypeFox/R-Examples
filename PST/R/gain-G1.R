## G1 gain function

G1 <- function(p1, p2, C, N, value=FALSE) {
	if (length(p1)!=length(p2)) {
		stop("p1 and p2 are not of equal lengths")
	} else if (any(p1==0) | any(p2==0)) {
		stop(" [!] one or more probability equals 0")
	}
	
	if (any(p1>=C*p2) || any(p1<=(1/C)*p2)) { gain <- TRUE } else { gain <- FALSE }

	return(gain)
}
