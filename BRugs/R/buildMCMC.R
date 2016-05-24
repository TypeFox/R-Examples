buildMCMC <- function(node, beg = samplesGetBeg(), end = samplesGetEnd(),
    firstChain = samplesGetFirstChain(), lastChain = samplesGetLastChain(), 
    thin = samplesGetThin()){
 
    oldBeg <- samplesGetBeg()
    oldEnd <- samplesGetEnd()
    oldFirstChain <- samplesGetFirstChain()
    oldLastChain <- samplesGetLastChain()
    oldThin <- samplesGetThin()
    on.exit({
        samplesSetBeg(oldBeg)
        samplesSetEnd(oldEnd)
        samplesSetFirstChain(oldFirstChain)
        samplesSetLastChain(oldLastChain)
        samplesSetThin(oldThin)
    })
    samplesSetBeg(beg)
    samplesSetEnd(end)
    samplesSetFirstChain(firstChain)
    samplesSetLastChain(lastChain)
    thin <- max(c(thin, 1))
    samplesSetThin(thin)
    mons <- samplesMonitors(node)
    
    subBuildMCMC <- function(node){
        sM <- samplesMonitors(node)
        if(length(sM) > 1 || sM != node)
            stop("node must be a scalar variable from the model, for arrays use samplesAutoC")
        sample <- samplesSample(node)
        numChains <- samplesGetLastChain() - samplesGetFirstChain() + 1
        matrix(sample, ncol = numChains)
    }

    sampleSize <- samplesSize(mons[1])
    end <- min(c(modelIteration(), samplesGetEnd()))
    thin <- samplesGetThin()
    numChains <- samplesGetLastChain() - samplesGetFirstChain() + 1
    sampleSize <- sampleSize %/% numChains
    beg <- end - sampleSize * thin + 1
    if (sampleSize==0) {
        mcmcobj <- NA
    }
    else { 
        samples <- lapply(mons, subBuildMCMC)
        samplesChain <- vector(mode="list", length=numChains)
        for(i in 1:numChains){
            if (is.R())
                temp <- sapply(samples, function(x) x[,i])
            else 
                temp <- sapply(samples, function(x,j) { x[,j]}, j=i)
##### If we want to special-case 1D-mcmc objects:
                                        #        if(ncol(temp) == 1){
                                        #            dim(temp) <- NULL
                                        #            samplesChain[[i]] <- temp
                                        #        }
                                        #        else{
            samplesChain[[i]] <- temp
            colnames(samplesChain[[i]]) <- mons
                                        #        }
        }
        mcmcobj <- lapply(samplesChain, mcmc, start = beg, end = end, thin = thin)
    }
    if(is.R())
        class(mcmcobj) <- "mcmc.list"
    else
        oldClass(mcmcobj) <- "mcmc.list"
    mcmcobj
}
