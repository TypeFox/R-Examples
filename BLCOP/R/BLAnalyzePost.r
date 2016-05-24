
###############################################################################
# Mango Solutions, Chippenham SN14 0SQ 2008
# posteriorFeasibility
# Author: Francisco
###############################################################################
# DESCRIPTION: Tries to assess the "feasibility" of a set of Black-Litterman views using the method described by Meucci and Fusai 
# in "Assessing Views".  This method is based on the Mahalanobis distance between the posterior and prior mean
# KEYWORDS: math
# TODO: Appears not to be completely correct at the moment
###############################################################################


posteriorFeasibility <- function(
    result                     # BLResult class object
)
{
    views <- result@views
    qv <- views@qv
    P <- views@P
    numAssets <- length(assetSet(views))
    sigmaInv <- solve(result@priorCovar)
    
    # calculates the Mahalanobis distance as described by the papaer
    mahal <- mahalanobis(result@posteriorMean, result@priorMean, cov = result@priorCovar, inverted = FALSE)
    mahalProb <- 1 - pchisq(mahal, df = numAssets)
    
    if(! result@kappa == 0)
      omega <-  diag(1 / views@confidences)
    else    
      omega <- result@kappa * tcrossprod(P %*% result@priorCovar, P)

    # 
    if(result@tau != 1)    
        warning("This function is not yet implemented for tau != 1, so the calculation of view senstivities will proceed assuming tau = 1")
   #  sensitivities <- -2 * dchisq(mahal, df = numAssets) * (solve(tcrossprod(P %*% result@priorCovar, P) 
   #     + omega) %*%  P %*%(result@posteriorMean-result@priorMean))    
    list("mahalDist" = mahal, "mahalDistProb" = mahalProb, sensitivities = "Not implemented yet") 
}