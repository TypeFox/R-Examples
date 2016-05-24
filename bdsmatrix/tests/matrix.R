library(bdsmatrix)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

tmat <- bdsmatrix(c(3,2,2,4), 
	      c(22,1,2,21,3,20,19,4,18,17,5,16,15,6,7, 8,14,9,10,13,11,12),
	      matrix(c(1,0,1,1,0,0,1,1,0,1,0,10,0,
                       0,1,1,0,1,1,0,1,1,0,1,0,10), ncol=2))
dimnames(tmat) <- list(NULL, letters[1:13])

smat <- as.matrix(tmat)
yy <- c(30,35,42,56,34,45,32,37,78,56,40,52,39)

# matrix multiplication
zz <- runif(13)
aeq(zz%*% smat, zz%*% tmat)
aeq(smat%*%zz, tmat%*% zz)

xx <- matrix(1:39, ncol=3)
aeq(smat %*% zz, tmat %*% zz)
aeq(t(xx) %*% smat, t(xx) %*% tmat)


amat <- tmat
amat@offdiag <- pi
bmat <- as.matrix(amat)

aeq(zz%*% amat, zz%*% bmat)
aeq(amat%*%zz, bmat%*% zz)


# Solve the right-hand side wrt a matrix
yy2 <- cbind(yy, -yy, yy+3)
zz1 <- solve(smat, yy2)
zz2 <- solve(tmat, yy2)
aeq(zz1, zz2)
aeq(zz2[,1], solve(tmat, yy))
