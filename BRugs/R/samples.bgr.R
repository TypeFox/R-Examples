"samplesBgr" <-
function(node, beg = samplesGetBeg(), end = samplesGetEnd(), 
    firstChain = samplesGetFirstChain(), lastChain = samplesGetLastChain(), 
    thin = samplesGetThin(), bins = 50, plot = TRUE, mfrow = c(3, 2), 
    ask = NULL, ann = TRUE, ...)
#   Plot bgr statistic                  
{
    mons <- samplesMonitors(node)
    if (any(grep("^inference can not be made", mons))) { stop(mons) }    
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
    beg <- max(beg, modelAdaptivePhase())
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
    result <- lapply(mons, plotBgr, bins = bins, plot = plot, ...)
    names(result) <- mons
    if(plot) invisible(result)
    else     return(result)
}
