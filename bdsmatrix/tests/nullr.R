#
# Make sure things work in the simpler case of no rmat component
#
library(bdsmatrix)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

tmat <- bdsmatrix(blocksize= c(3,2,2,4), 
              blocks= c(22,1,2,21,3,20,19,4,18,17,5,16,15,6,7, 
                        8,14,9,10,13,11,12))
dimnames(tmat) <- list(NULL, letters[1:11])
smat <- as.matrix(tmat)
yy <- c(30,35,42,56,34,45,32,37,78,56,40)

aeq(diag(tmat), diag(smat))
zz <- seq(1,11,2)
aeq(as.matrix(tmat[zz,zz]), smat[zz,zz])

ch0 <- chol(smat)
ch1 <- gchol(smat)
ch2 <- gchol(tmat)

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
