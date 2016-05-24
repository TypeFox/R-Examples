para_trans <-
function(transmat) {
M = dim(transmat)[1]
A = matrix(0,M,M-1)
for ( i in 1:(M-1)) {
    A[,i] = log(transmat[,i]/transmat[,M])
	}
	return(A)
}
