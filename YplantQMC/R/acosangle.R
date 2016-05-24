
# From Yplant  (leaf.makeviewbound in Dpunit.pas)
acosangle <- function (a, b) {
    alen <- sqrt(sum(a^2))
    blen <- sqrt(sum(b^2))
	
	acosnum <- abs(sum(a*b))
	
	if(alen > 0)ang <- acosnum / alen
	if(blen > 0)ang <- acosnum / blen
	
    return(ang)
}
