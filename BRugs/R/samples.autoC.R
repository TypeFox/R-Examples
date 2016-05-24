"samplesAutoC" <-
function(node, chain, beg = samplesGetBeg(), end = samplesGetEnd(),
thin = samplesGetThin(), plot = TRUE, mfrow = c(3, 2), ask = NULL, ann = TRUE, ...)
#   Plot auto correlation function
{
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
    chain <- max(c(1, chain))
    chain <- min(c(getNumChains(), chain))
    samplesSetFirstChain(chain)
    samplesSetLastChain(chain)
    thin <- max(c(thin, 1))
    samplesSetThin(thin)
    mons <- samplesMonitors(node)
    if(plot){
        if (is.R())
        par(mfrow = mfrow, ask = ask, ann = ann)
        else
        par(mfrow = mfrow, ask = ask)
    }
    result <- lapply(mons, plotAutoC, plot = plot, ...)
    names(result) <- mons
    if(plot) invisible(result)
    else     return(result)
}
