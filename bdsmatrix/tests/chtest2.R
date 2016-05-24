#
# Inverse of the matrix:  
#
library(bdsmatrix)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

tmat <- bdsmatrix(c(3,2,2,4), 
	      c(22,1,2,21,3,20,19,4,18,17,5,16,15,6,7, 8,14,9,10,13,11,12),
	      matrix(c(1,0,1,1,0,0,1,1,0,1,0,10,0,
                       0,1,1,0,1,1,0,1,1,0,1,0,10), ncol=2))
dimnames(tmat) <- list(NULL, letters[1:13])

smat <- as.matrix(tmat)

inv1 <- solve(smat)
inv2 <- as.matrix(solve(tmat)) # the result is a full, non-sparse matrix
aeq(inv1, inv2)

inv3 <- solve(gchol(tmat))         #sparse version, not all parts will be there
inherits(inv3, 'bdsmatrix')         #This should be true
aeq(inv3@blocksize, tmat@blocksize) # Should be the same shape at tmat
inv3 <- as.matrix(inv3)		    # What is returned should be correct
aeq(inv1[1:3,1:3], inv3[1:3, 1:3])
aeq(inv1[4:5,4:5], inv3[4:5, 4:5])
aeq(inv1[6:7,6:7], inv3[6:7, 6:7])
aeq(inv1[8:11,8:11], inv3[8:11, 8:11])
aeq(inv1[,12:13], inv3[, 12:13])    # and rmat the same too
