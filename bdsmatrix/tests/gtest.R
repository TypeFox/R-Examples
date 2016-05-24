library(bdsmatrix)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

tmat <- bdsmatrix(c(3,2,2,4), 
	      c(22,1,2,21,3,20,19,4,18,17,5,16,15,6,7, 8,14,9,10,13,11,12),
	      matrix(c(1,0,1,1,0,0,1,1,0,1,0,10,0,
                       0,1,1,0,1,1,0,1,1,0,1,0,10), ncol=2))
dimnames(tmat) <- list(NULL, letters[1:13])

smat <- as.matrix(tmat)

# Create a matrix that is symmetric, but not positive definite
#   The first one, temp, has column 6 redundant with cols 1-5
temp <- smat[c(1:5, 5:10), c(1:5, 5:10)]
ch1  <- gchol(temp)
aeq(diag(ch1)[6], 0)  # Check that it has a zero in the proper place
ginv <- solve(ch1)    # see if I get a generalized inverse
aeq(temp %*% ginv %*% temp, temp)
aeq(ginv %*% temp %*% ginv, ginv)

# Now create one that is negative definite 
ch2 <- gchol(smat)
temp2 <- as.matrix(ch2)
temp3 <- diag(ch2) * rep(c(1, -1), length=nrow(smat))
xmat  <- temp2 %*% diag(temp3) %*% t(temp2)
xmat  <- (xmat + t(xmat))/2  #work out round-off errors
ch3 <- gchol(xmat)

aeq(diag(ch3), temp3)
aeq(as.matrix(ch3), temp2)
