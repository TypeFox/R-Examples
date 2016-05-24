library(pcalg)

set.seed(34)
g <- randomDAG(5,0.8)
trMat <- matrix(0, 5,5)
trMat[1,5] <- 0.305
trMat[1,4] <- 0.863
trMat[1,2] <- 0.354
trMat[2,4] <- 0.392
trMat[2,5] <- 0.495
trMat[3,4] <- 0.278


## eff 1->5: 0.305 + 0.354*0.495
trEff <- 0.305 + 0.354*0.495
estEff <- causalEffect(g, 5,1)

if (!(round(trEff,3) == round(estEff, 3))) {
  stop("Test in wgtMatrix: Weights have wrong value!")
}
