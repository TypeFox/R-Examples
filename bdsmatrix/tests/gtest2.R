library(bdsmatrix)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

#
# Test multiplication of a vector/matrix times a gchol
#
tmat <- bdsmatrix(c(3,2,2,4), 
              c(22,1,2,21,3,20,19,4,18,17,5,16,15,6,7, 8,14,9,10,13,11,12),
              matrix(c(1,0,1,1,0,0,1,1,0,1,0,10,0,
                       0,1,1,0,1,1,0,1,1,0,1,0,10), ncol=2))
dimnames(tmat) <- list(NULL, letters[1:13])

gmat <- gchol(tmat)
g2 <- as.matrix(gmat) %*% diag(sqrt(diag(gmat)))


aeq(1:13 %*% g2, 1:13 %*% gmat)  #vectors first
aeq(g2 %*% 1:13, gmat %*% 1:13)

temp <- matrix(runif(39), nrow=3)
aeq(temp %*% g2, temp %*% gmat)
aeq(g2 %*% t(temp), gmat %*% t(temp))
