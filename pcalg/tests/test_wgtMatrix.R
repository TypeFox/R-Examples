library(pcalg)

set.seed(42)

## wmat_ij is edge from j to i
g <- randomDAG(3,0.4)
wmat <- wgtMatrix(g)
if (!(wmat[2,1]!=0 & wmat[1,2]==0 & wmat[3,1]!=0 & wmat[3,2]==0)) {
  stop("Test of wgtMatrix: Something with orientation of edges is wrong!")
}

## test weird parameters
g <- randomDAG(3,0)
wmat <- wgtMatrix(g)
if (!all(wmat==matrix(0,3,3))) {
  stop("Test of wgtMatrix: Problem when used on empty graph!")
}

## test if weight parameters are correct
set.seed(34)
g <- randomDAG(5,0.8)
trMat <- matrix(0, 5,5)
trMat[1,5] <- 0.305
trMat[1,4] <- 0.863
trMat[1,2] <- 0.354
trMat[2,4] <- 0.392
trMat[2,5] <- 0.495
trMat[3,4] <- 0.278

if (!all(round(wgtMatrix(g),3) == t(trMat))) {
  stop("Test in wgtMatrix: Weights have wrong value!")
}
