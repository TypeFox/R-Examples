#
# Test out the behavior of a 1x1 bds corner.
#  Actually, the problem that motivated this occurred whenever the
# rmat portion was larger than the block diagonal portion.
#
library(bdsmatrix)
test1 <- bdsmatrix(blocksize=1, blocks=33, 
                   rmat=matrix(c(17,33,7,-1, -7,7,48,-7,
                                 1, -1, -7,4),4))

test2 <- bdsmatrix(blocksize=2, blocks=c(33,17,33),
                    rmat=matrix(c( -7,7,48,-7, 1, -1, -7,4),4))
all.equal(as.matrix(test1), as.matrix(test2))

g1 <- gchol(test1)
g2 <- gchol(test2)
all.equal(as.matrix(g1), as.matrix(g2))

s1 <- solve(g1, full=T)
s2 <- solve(g2, full=T)
all.equal(as.matrix(s1), as.matrix(s2))

all.equal(solve(test1), solve(test2))
