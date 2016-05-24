library(bdsmatrix)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

tmat <- bdsmatrix(c(3,2,2,4), 
	      c(22,1,2,21,3,20,19,4,18,17,5,16,15,6,7, 8,14,9,10,13,11,12),
	      matrix(c(1,0,1,1,0,0,1,1,0,1,0,10,0,
                       0,1,1,0,1,1,0,1,1,0,1,0,10), ncol=2))
dimnames(tmat) <- list(NULL, letters[1:13])

smat <- as.matrix(tmat)
yy <- c(30,35,42,56,34,45,32,37,78,56,40,52,39)
# look at inverses more closely
#  (I needed this when some of the other tests weren't being passed,
#  to figure out where in the decomposition/inversion/multiply process
#  the flaw was).

ch1 <- gchol(smat)
ch2 <- gchol(tmat)

inv1 <- solve(as.matrix(ch1))
inv2 <- solve(ch2,full=F)  #inverse of the cholesky, not of tmat
aeq(inv1, as.matrix(inv2))

# Full matrix tests
inv3 <- solve(smat)
inv4 <- solve(tmat)
inv5 <- solve(gchol(smat), full=T)
aeq(inv3, inv4)
aeq(inv3, inv5)

# The following test is false by design: when called with a bdsmatrix
#  object that has an rmat portion, the true inverse is dense.  But
#  coxme only needs the trace for one calcluation; solve(gchol(tmat))
#  cheats and only returns the block diagonal portion of the inverse.
#inv6 <- solve(gchol(tmat), full=T)
#aeq(inv3, inv6)

#
# Now test the solution to a partial solve
#  We want to be able to transform a matrix to uncorrelated form
#  If tmat= LDL', and A is general, I want (D^{-1/2}) L^{-1} A
#
amat <- matrix(runif(5*nrow(tmat)), nrow=nrow(tmat))
xx1 <- diag(1/sqrt(diag(ch1))) %*% solve(as.matrix(ch1), amat) 
xx2 <- solve(ch2, amat, full=F)
aeq(xx1, xx2)

xx1 <- diag(1/sqrt(diag(ch1))) %*% solve(as.matrix(ch1), yy)
xx2 <- solve(ch2, yy, full=F)
aeq(xx1, xx2)
