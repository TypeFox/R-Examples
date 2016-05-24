## G2 gain function

G2 <- function(p1, p2, C, N, value=FALSE) {
	if (length(p1)!=length(p2)) {
		stop("p1 and p2 are not of equal lengths")
	} else if (any(p1==0) | any(p2==0)) {
		stop(" [!] one or more probability equals 0")
	} 

	delta <- sum(p1*log(p1/p2))
	if ((delta*N) < C) { gain <- FALSE } else { gain <- TRUE }

	if (value) { return(delta) } else { return(gain) }
}
