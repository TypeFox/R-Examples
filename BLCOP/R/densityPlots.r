###############################################################################
# Mango Solutions, Chippenham SN14 0SQ 2008
# densityPlots
# Author: Francisco
###############################################################################
# DESCRIPTION: Generic function that displays density plots of the prior and posterior distributions
# KEYWORDS: hplot
###############################################################################


densityPlots <- function(result, assetsSel = NULL, numSimulations = BLCOPOptions("numSimulations"), ...)
{
    stop("Not implemented for this class")
}

setGeneric("densityPlots")

densityPlots.COPViews <- function(result, assetsSel = seq(along = result@views@assets) , numSimulations = BLCOPOptions("numSimulations") ,...) 
{
    marketSims <- sampleFrom(result@marketDist, numSimulations)[,drop=FALSE]
    colnames(marketSims) <- assetSet(result@views)
    for(i in seq(along = assetsSel))
    {
        sims <- tail(result@posteriorSims[, assetsSel[i]], numSimulations)
		plot(density(sims), col = "blue", xlab = result@views@assets[assetsSel[i]], 
            main = "Kernel density estimates of posterior and prior", ...)
        lines(density(marketSims[,assetsSel[i]]), col = "black", ...)
        abline(v = mean(sims, lty = 2, col = "blue"))
        abline(v = mean(marketSims[,i]), lty = 2, col = "black")
    }
    legend(x = "topright", legend = c("Posterior", "Prior"), lty = c(1,1), col = c("blue", "black"))
}

setMethod("densityPlots", signature(result = "COPResult"), densityPlots.COPViews)

densityPlots.BLResult <- function(result, assetsSel = seq(along = assetSet(result@views)) , numSimulations = BLCOPOptions("numSimulations"),
                                    ...)
{
    for(i in seq(along = assetsSel))
    {
        postMean <- result@posteriorMean[assetsSel[i]] 
        priorMean <- result@priorMean[assetsSel[i]]
        postStDev <- sqrt(result@posteriorCovar[assetsSel[i],assetsSel[i]] )
        priorStDev <- sqrt(result@priorCovar[assetsSel[i],assetsSel[i]])
      
        plotDispersion <- max(postStDev, priorStDev)
        x <- seq(from = min(priorMean,postMean) - 2.5 * abs(plotDispersion), to = max(priorMean,postMean) + 2.5 * abs(plotDispersion), length = 200)
        xLabel <- if(is.character(assetsSel)) assetsSel[i] else result@views@assets[i]
        
        if(dnorm(postMean, mean = postMean, sd = postStDev) < dnorm(priorMean, mean = priorMean, sd = priorStDev))
        {                
            plot(x, dnorm(x, mean = priorMean, sd = priorStDev), col = "black", type = "l",..., ylab = "Density", xlab = xLabel)
            abline(v = priorMean, lty = 2, col = "black")        
            lines(x, dnorm(x, mean = postMean, sd = postStDev), col = "blue", type = "l",...)
            abline(v = postMean, lty = 2, col = "blue")
            legend(x = "topright", legend = c("Prior", "Posterior"), lty = c(1,1), col = c("black", "blue"))
        }
        else
        {                
            plot(x, dnorm(x, mean = postMean, sd = postStDev), col = "blue", type = "l",..., ylab = "Density", xlab = xLabel)
            abline(v = postMean, lty = 2, col = "blue")

            lines(x, dnorm(x, mean = priorMean, sd = priorStDev), col = "black", type = "l",...)
            abline(v = priorMean, lty = 2, col = "black")        
            legend(x = "topright", legend = c("Prior", "Posterior"), lty = c(1,1), col = c("black", "blue"))   
        }
    }      
}

setMethod("densityPlots", signature(result = "BLResult"), densityPlots.BLResult)


biDensityPlots <- function(result, assetsSel , numSimulations = BLCOPOptions("numSimulations"), nBins,
		...)
{
	.assertClass(result, "COPResult")
	stopifnot(length(assetsSel) == 2)
	marketSims <- sampleFrom(result@marketDist, numSimulations) 
	colnames(marketSims) <- assetSet(result@views)
	
	marketSims <- marketSims[,assetsSel,drop=FALSE]
	
	sims <- tail(result@posteriorSims, numSimulations)[,assetsSel]
	
	hexBin <- hexBinning(sims, bins = nBins)
	
	par(mfrow = c(1,2))
	plot(hexBin,  xlab = assetsSel[1], ylab = assetsSel[2], main = "Posterior", col = rev(greyPalette(nBins)))
	
	hexBin <- hexBinning(marketSims, bins = nBins)
	plot(hexBin, xlab = assetsSel[1], ylab = assetsSel[2], main = "Prior", col = seqPalette(nBins))
	
}