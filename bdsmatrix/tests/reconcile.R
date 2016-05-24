#
# Test cases for bdsmatrix.reconcile
#
library(bdsmatrix)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

tmat <- bdsmatrix(c(3,2,2,4), 
              c(22,1,2,21,3,20,19,4,18,17,5,16,15,6,7, 8,14,9,10,13,11,12),
              matrix(c(1,0,1,1,0,0,1,1,0,1,0,10,0,
                       0,1,1,0,1,1,0,1,1,0,1,0,10), ncol=2))
id <- letters[1:13]
dimnames(tmat) <- list(id, id)

rmat <- matrix(1:169, 13,13)
rmat <- (rmat + t(rmat))/2
dimnames(rmat) <- list(rev(id), rev(id))

xmat <- bdsmatrix.reconcile(list(tmat, rmat), group=id)
temp <- xmat[[1]]
aeq(temp@blocksize, 13)
aeq(as.matrix(temp), as.matrix(tmat))

temp <- xmat[[2]]
aeq(temp@blocksize, 13)
aeq(as.matrix(temp), rmat[13:1, 13:1])

xmat <- bdsmatrix.reconcile(list(rmat, bdsI, tmat), group= rev(id))
temp <- xmat[[1]]
aeq(as.matrix(temp), rmat)

temp <- xmat[[2]]
aeq(as.matrix(temp), diag(13))

temp <- xmat[[3]]
aeq(as.matrix(temp), (as.matrix(tmat))[13:1,13:1])

#
# Simplest case
#
xmat <- bdsmatrix.reconcile(bdsI, id)
all(xmat@blocksize==1)
aeq(as.matrix(xmat), diag(13))

#
# The case list(tmat, tmat) will fail  - can't have 2 rmats
#
#xmat <- bdsmatrix.reconcile(list(tmat, tmat), id)

xmat <- bdsmatrix.reconcile(list(tmat, bdsI), rev(id))
temp <- xmat[[1]]
aeq(as.matrix(temp), as.matrix(tmat))
aeq(dimnames(temp)[[1]], id)

aeq(as.matrix(xmat[[2]]), diag(13))

#
# Now for the hard one: 2 bdsmatrices, different orders, different
#   blocksize, but one contains the other
#
tmat <- tmat[1:11, 1:11]
tord <- c(11:8, 2,1,3, 6,7,5,4)
rmat <- (as.matrix(tmat))[tord, tord]
rmat <- bdsmatrix(blocksize=c(4,5,2),
                  blocks=c(rmat[1:4,1:4], rmat[5:9,5:9], rmat[10:11, 10:11]),
                  dimnames=list(id[tord], id[tord]))

aeq(as.matrix(tmat)[tord,tord], as.matrix(rmat))

xmat <- bdsmatrix.reconcile(list(tmat/2, rmat), id[1:11])
all.equal(xmat[[2]], rmat)
all.equal(xmat[[1]]*2, rmat)

# Now toss out a row/col
#  Give it a different name, too
xx <- id[tord]
xx[1] <- 44
dimnames(rmat) <- list(xx,xx)
xmat <- bdsmatrix.reconcile(list(tmat/2, rmat), id[1:9])
all.equal(xmat[[1]]*2, xmat[[2]])
