#The following code implements the calculations of ENMtools and was
# created by Hackett in the following google groups mailing
# list: https://groups.google.com/forum/#!topic/maxent/EsXZKlpvdTI

MaxentIC <- function(csvfile, grdfile, lambdasfile) {
    nparams = 0
    probsum = 0
    loglikelihood = 0
    AICcscore = 0
    AICscore = 0
    BICscore = 0
    
    lambdases <- read.csv(lambdasfile, header=FALSE)
    nparams <- nrow(lambdases[lambdases$V2 != 0, ])
    nparams = nparams - 4
    
    layerRaw <- raster::raster(grdfile)
    probsum <- raster::cellStats(layerRaw, sum)
    
    points <- read.csv(csvfile)
    npoints <- nrow(points)
    layerValues <- raster::extract(layerRaw, points[, c("longitude", "latitude")])
    loglikelihood <- sum(log(layerValues / probsum))
    
    if (nparams >= npoints - 1) {
        AICcscore <- "x"
        AICscore <- "x"
        BICscore <- "x"
    } else {
        AICcscore = (2 * nparams - 2 * loglikelihood) + (2 * (nparams) * (nparams + 1) / (npoints - nparams - 1))
        AICscore = 2 * nparams - 2 * loglikelihood
        BICscore = nparams * log(npoints) - 2 * loglikelihood
    }
    
    ICs <- c(npoints, nparams,  loglikelihood, AICscore, AICcscore, BICscore)
    
    return(ICs)
}
