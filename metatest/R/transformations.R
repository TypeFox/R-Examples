# fisher's z-transformation
r2z <-
function(r) {
	return(0.5*(log(1+r)-log(1-r)))
}

z2r <-
function(z) {
	return( (exp(2*z)-1) / (exp(2*z)+1) )
}

z2d <-
function(z) {
	return(r2d(z2r(z)))
}

r2d <-
function(r) {
	return(2*r/(sqrt(1-r*r)))
}

