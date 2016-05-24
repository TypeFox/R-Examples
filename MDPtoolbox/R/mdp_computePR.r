# function PR = mdp_computePR(P,R)

mdp_computePR <- function(P,R) {

# Test if R has form R(SxA)
if (length(dim(R)) == 2) { # & ~iscell(R)
	PR <- R
} else { # else R has form R(SxSxA)
	if (is.list(R)){
		S<-dim(R[[1]])[1]
	}else{
		S <- dim(R)[1]
	}
	if (is.list(P)) {
		A <- length(P)
		PR <- array(NA, c(S,A))
		if (is.list(R)) {
			for (a in 1:A) {
				PR[,a] <- rowSums(P[a][[1]]*R[a][[1]])
			}
		} else {
			for (a in 1:A) {
				PR[,a] <- rowSums(P[a][[1]]*R[,,a])
			}
		}
	} else {
		A <- dim(P)[3]
		PR <- array(NA, c(S,A))
		if (is.list(R)) {
			for (a in 1:A) {
				PR[,a] <- rowSums(P[,,a]*R[a][[1]])
			}
		} else {
			for (a in 1:A) {
				PR[,a] <- rowSums(P[,,a]*R[,,a])
			}
		}
	}
}

return(PR)

}
