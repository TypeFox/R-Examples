###############################################################################
# Mango Solutions, Chippenham SN14 0SQ 2008
# COPPosterior
# Author: Francisco
###############################################################################
# DESCRIPTION: Calculates the posterior distribution of "the market" given a set of views.  The posterior is returned
# in the form of a set of simulations.
# KEYWORDS: math
###############################################################################


COPPosterior <- function
(
    marketDist,                                       # mvdistribution object that defines the distribution of "the market"
    views,                                            # COPViews object
    numSimulations = BLCOPOptions("numSimulations")  # number of samples to use for each monte-carlo simulation
) 
{  
    
    # generate simulations of the market distribution and order them. 
    marketSimulations <- t(sampleFrom(marketDist, numSimulations))
    
    # now simulate from each subjective view
    subjSimulations <- sapply(views@viewDist, sampleFrom, n = numSimulations )
    numViews <- length(views@viewDist)
    
    # compute the orthogonal complement of the pick matrix
    nullPick <- t(Null(t(views@pick)))
    pick <- views@pick
    # calculate the product of "the market" with the pick matrix 
    impliedViews <- pick %*% marketSimulations
    # calculate the orthongal complement of the above productg
    complement <- nullPick %*% marketSimulations
    
    # Now generate samples from blended views and "implied market views".
    
    .innerChoiceSample <- function(conf) {
        sample(0:1, prob = c(1-conf, conf), numSimulations, replace = TRUE)
    }
    choices <- t(sapply(views@confidences, .innerChoiceSample))
    combinedSimulations <- matrix(0, nrow = numViews, ncol = numSimulations)
    combinedSimulations[choices == 0] <- impliedViews[choices == 0]
    combinedSimulations[choices == 1] <- t(subjSimulations)[choices==1]
    # combinedSimulations <- (1-views@confidences) * impliedViews + views@confidences * t(subjSimulations)

    impliedCopula <- array(dim = dim(impliedViews))
    pooledSimulations <- array(dim = dim(combinedSimulations))
    # compute the copula of the implied views
    for(i in 1:nrow(impliedViews))
    {
        cdf <- .empCDF(impliedViews[i,])
        impliedCopula[i,] <- cdf(impliedViews[i,])
        quant <- .empQuantile(combinedSimulations[i,])
        pooledSimulations[i,] <- quant(impliedCopula[i,])       
    }
    # rotate back to "market coordinates"
    rotMatrix <- solve(rbind(pick, nullPick))
    result <- t(rotMatrix %*% rbind(pooledSimulations, complement))
    colnames(result) <- assetSet(views)
    new("COPResult", views = views, marketDist = marketDist, posteriorSims = result)
    
}
                 