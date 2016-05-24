library(HyperbolicDist)
# pi, zeta, delta, mu
hypDat <- rhyperb(1000, c(0, 1, 1, 0))
hResultList <- hyperbFitStart(hypDat)
hResult <- hResultList$ThetaStart
hResult <- c(hResult[1], exp(hResult[2]), exp(hResult[3]), hResult[4])
hResult <- hyperbChangePars(1, 2, hResult)
hResult <- c(hResult[4], hResult[3],
             hResult[1], hResult[2])
names(hResult) <- c("mu", "delta", "alpha", "beta")
hResultList$ThetaStart <- hResult
hResult <- hResultList
detach("package:HyperbolicDist")
library(GeneralizedHyperbolic)
# should be close to c(0, 1, 0, 1)
gResult <- hyperbFitStart(hypDat)
detach("package:GeneralizedHyperbolic")


cat("Test results:\n\n", "HyperbolicDist result:\n")
print(hResult$ThetaStart)
cat("GeneralizedHyperbolic result:\n")
print(gResult$paramStart)
