#
# Test out conversions between bdsmatrix and Matrix objects
#
library(coxme)
require(bdsmatrix)
require(Matrix)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

tmat <- bdsmatrix(c(3,2,2,4), 
	      c(22,1,2,21,3,20,19,4,18,17,5,16,15,6,7, 8,14,9,10,13,11,12),
	      matrix(c(1,0,1,1,0,0,1,1,0,1,0,10,0,
                       0,1,1,0,1,1,0,1,1,0,1,0,10), ncol=2))
dimnames(tmat) <- list(NULL, letters[1:13])

smat <- as.matrix(tmat)
cmat <- as(tmat, 'dsCMatrix')
cmat2 <- as(forceSymmetric(smat), 'dsCMatrix')

aeq(as.matrix(smat), as.matrix(cmat))
aeq(as.matrix(tmat), as.matrix(cmat2))

tmat2 <- as(cmat, 'bdsmatrix')
aeq(smat, as.matrix(tmat2))

# The above makes a 13 by 13, not a very elegant nor complete test
#  This one is more strict
# This fails with Matrix 3.1.3 due to an error in the Matrix package -
#  subscripting can add row names.  Add it back at a later date.
if (FALSE) {
    tmat3 <- as(cmat[1:11, 1:11], 'bdsmatrix')
    all.equal(tmat3, tmat[1:11, 1:11])
}
