#
# Test out the Cholesky, special case of a null block corner 
# In this case there is no advantage to a bdsmatrix as it consists of only
#  the ordinary matrix part.  This case arises in coxme with an (x1+x2 | 1)
#  term, however, so it is nice to have it work instead of coding lots of
#  if/else logic in that code base.
#
library(bdsmatrix)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

tmat <- bdsmatrix(c(3,2,2,4), 
	      c(22,1,2,21,3,20,19,4,18,17,5,16,15,6,7, 8,14,9,10,13,11,12),
	      matrix(c(1,0,1,1,0,0,1,1,0,1,0,10,0,
                       0,1,1,0,1,1,0,1,1,0,1,0,10), ncol=2))
dimnames(tmat) <- list(NULL, letters[1:13])
smat <- as.matrix(tmat)
tmat <- bdsmatrix(integer(0), numeric(0), rmat=smat)  
yy <- c(30,35,42,56,34,45,32,37,78,56,40,52,39)

aeq(diag(tmat), diag(smat))
zz <- seq(1,13,2)
aeq(as.matrix(tmat[zz,zz]), smat[zz,zz])

ch0 <- chol(smat)
ch1 <- gchol(smat)
ch2 <- gchol(tmat)
# The gchol routines use the composition LDL', where L is lower triangular
#  with a diagonal of 1's, and D is diagonal.  chol() uses U'U where U is
#  upper trangular.  
# The as.matrix function returns L and the diag function returns D.
#  Convert and compare
aeq(diag(ch1), diag(ch2))
temp <- as.matrix(ch2)
aeq(temp, as.matrix(ch1))
temp3 <- temp %*% diag(sqrt(diag(ch2))) 
aeq(temp3, t(ch0))

zz0 <- solve(smat, yy)
zz1 <- solve(ch1, yy)
zz2 <- solve(tmat, yy)
aeq(zz1, zz2)
aeq(zz0, zz1)

inv1 <- solve(smat)
inv2 <- as.matrix(solve(tmat)) # the result is a full, non-sparse matrix
aeq(inv1, inv2)

inv3 <- solve(gchol(tmat))
aeq(inv1, as.matrix(inv3))

gmat <- gchol(tmat)
g2 <- as.matrix(gmat) %*% diag(sqrt(diag(gmat)))
aeq(1:13 %*% g2, 1:13 %*% gmat)  #vectors first
aeq(g2 %*% 1:13, gmat %*% 1:13)
temp <- matrix(runif(39), nrow=3)
aeq(temp %*% g2, temp %*% gmat)
aeq(g2 %*% t(temp), gmat %*% t(temp))
