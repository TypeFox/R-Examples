ojaMedianControl <- function(sigmaInit = 0.0, sigmaAda = 20, adaFactor = 0.5, iter = 1e+06, useAllSubsets = FALSE, nSubsetsUsed = 1000,
                              sigmaLog10Dec = 10, storeSubDet = TRUE, eps = 0.1, chi2 = 0.95, 
                              samples = 20, maxlines=1000, S1 = cov, S2 = cov4, S1args = list(), S2args = list(), volume = 1e-6, boundedExact = T)
{
list(sigmaInit=sigmaInit, sigmaAda=sigmaAda, adaFactor=adaFactor,iter=iter,useAllSubsets=useAllSubsets, nSubsetsUsed=nSubsetsUsed,
     sigmaLog10Dec=sigmaLog10Dec, storeSubDet=storeSubDet, eps=eps, chi2=chi2, 
     samples=samples, maxlines=maxlines, S1 = S1, S2 = S2, S1args = S1args, S2args = S2args, volume = volume, boundedExact = boundedExact)
}
