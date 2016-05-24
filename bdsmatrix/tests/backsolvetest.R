library(bdsmatrix)
#
# A test of the backsolve function
#
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

tmat <- matrix(rep(1:5,5), 5, 5)
tmat <- tmat + t(tmat)
diag(tmat) <- diag(tmat) + 10

gt <- gchol(tmat)
g1 <- as.matrix(gt)
gd <- diag(sqrt(diag(gt)))
gc <- gd %*% t(g1)  #usual cholesky form

xmat <- cbind(1:5, 11:15)

s1 <- backsolve(gt, xmat, upper=TRUE)  #the default
aeq(gd %*% t(g1) %*% s1, xmat)
all.equal(s1, backsolve(gc, xmat)) 

s2 <- backsolve(gt, xmat, upper=FALSE)
aeq(g1 %*% gd %*% s2, xmat)
all.equal(backsolve(gt,xmat, upper=F), backsolve(t(gc),xmat, upper=F)) 


# Now for bdsmatrix objects
tmat <- bdsmatrix(c(3,2,2,4), 
	      c(22,1,2,21,3,20,19,4,18,17,5,16,15,6,7, 8,14,9,10,13,11,12),
	      matrix(c(1,0,1,1,0,0,1,1,0,1,0,10,0,
                       0,1,1,0,1,1,0,1,1,0,1,0,10), ncol=2))
dimnames(tmat) <- list(NULL, letters[1:13])
smat <- as.matrix(tmat)

gt <- gchol(tmat)
gs <- gchol(smat)

xmat <- cbind(1:13, 1:13*2 + 3)

s1 <- backsolve(gt, xmat)
s2 <- backsolve(gs, xmat)
s3 <- backsolve(gt, xmat, upper=FALSE)
s4 <- backsolve(gs, xmat, upper=FALSE)

aeq(s1, s2)
aeq(s3, s4)

g1 <- as.matrix(gt)
gd <- diag(sqrt(diag(gt)))
aeq(gd %*% t(g1) %*% s1, xmat)
aeq(g1 %*% gd %*% s3, xmat)
