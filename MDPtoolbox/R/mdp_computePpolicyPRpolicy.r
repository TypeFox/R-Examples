mdp_computePpolicyPRpolicy <- function(P, R, policy) {

if (is.list(P)) {
	S <- dim(P[[1]])[1]
} else {
	S <- dim(P)[1]
}
Ppolicy <- matrix(0,S,S)
PRpolicy <- numeric(S)

if (is.list(P)) {
	A <- length(P)
} else {
	A <- dim(P)[3]
}

for (a in 1:A) {
	ind <- which(policy == a)
	if (length(ind) > 0) {
		if (is.list(P)) {
			Ppolicy[ind,] <- as.matrix(P[[a]][ind,])
		} else {
			Ppolicy[ind,] <- P[ind,,a]
		}
		PR <- mdp_computePR(P,R)
		PRpolicy[ind] <- PR[ind,a]
# Matlab said PRpolicy(ind,1) = PR(ind,a); but this seems to be OK for a vector
	}
}

if ( inherits(PR, 'sparseMatrix-class') ) {
	PRpolicy <- Matrix(PRpolicy, sparse = T)
}

return(list(Ppolicy, PRpolicy))

}
