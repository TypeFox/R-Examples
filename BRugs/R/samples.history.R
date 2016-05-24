"samplesHistory" <-
function(node, beg = samplesGetBeg(), end = samplesGetEnd(),
firstChain = samplesGetFirstChain(), lastChain = samplesGetLastChain(), 
thin = samplesGetThin(), plot = TRUE, mfrow = c(3, 1), ask = NULL, ann = TRUE, ...)
#   Plot history                    
{
    sM <- samplesMonitors(node)[1]
    if(sM == "model must be initialized before monitors used")
        stop("model must be initialized / updated / monitored before samplesSample is used")
    if(length(grep("^no monitor set for variable", sM)))
        stop(sM)
    if (samplesSize(sM[1])==0) 
        stop("No monitored samples available")
    if(plot && is.null(ask)) {
      if (is.R())
        ask <- !((dev.cur() > 1) && !dev.interactive())
      else
        ask <- !((dev.cur() > 1) && !interactive())
    }
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
    if(plot){
        if (is.R())
        par(mfrow = mfrow, ask = ask, ann = ann)
        else 
        par(mfrow = mfrow, ask = ask)
    }
    result <- lapply(mons, plotHistory, plot = plot, ...)
    names(result) <- mons
    if(plot) invisible(result)
    else     return(result)

}
