#
# Test out math aspects
#
library(bdsmatrix)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

tmat <- bdsmatrix(c(3,2,2,4), 
	      c(22,1,2,21,3,20,19,4,18,17,5,16,15,6,7, 8,14,9,10,13,11,12),
	      matrix(c(1,0,1,1,0,0,1,1,0,1,0,10,0,
                       0,1,1,0,1,1,0,1,1,0,1,0,10), ncol=2))
dimnames(tmat) <- list(NULL, letters[1:13])

smat <- as.matrix(tmat)

aeq(smat+2.1, as.matrix(tmat+2.1))
aeq(smat/2.1, as.matrix(tmat/2.1))
aeq(smat-2.1, as.matrix(tmat-2.1))
aeq(smat*2.1, as.matrix(tmat*2.1))
aeq(round(smat,1), as.matrix(round(tmat,1)))
aeq(exp(smat), as.matrix(exp(tmat)))

aeq(sum(smat), sum(tmat))
aeq(prod(smat), prod(tmat))
aeq(sum(smat+3), sum(tmat+3))
aeq(prod(smat+2), prod(tmat+2))
aeq(range(smat), range(tmat))
aeq(max(smat), max(tmat))
aeq(min(smat), min(tmat))

aeq(smat+1:13, tmat+1:13)
aeq(smat+1:13, 1:13 +tmat )
aeq(smat+tmat, 2*smat)
all.equal(tmat+tmat, 2*tmat)

aeq(sort(unique(c(smat))), sort(unique(tmat)))

#
# check out the alternate input style, with full blocks
#
rmat <- bdsmatrix(c(3,2,2,4),
		  c(22,1,2,1,21,3,2,3,20, 19,4,4,18, 17,5,5,16,
		    15,6,7,8,6,14,9,10,7,9,13,11,8,10,11,12),
		  matrix(c(1,0,1,1,0,0,1,1,0,1,0,10,0,
			   0,1,1,0,1,1,0,1,1,0,1,0,10), ncol=2),
		  dimnames=list(NULL, letters[1:13]))
all.equal(rmat, tmat)


# Do some subscripting
zz <- c(1,2,7,8,9,11)
aeq(smat[zz,zz], as.matrix(tmat[zz,zz]))

all.equal(smat[zz, 8:13], tmat[zz, 8:13])  # both are matrices

# Diagonals
aeq(diag(smat), diag(tmat))
zz <- diag(smat)
diag(smat) <- zz*2
diag(tmat) <- zz*2
all.equal(smat, as.matrix(tmat))
