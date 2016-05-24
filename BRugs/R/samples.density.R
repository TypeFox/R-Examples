"samplesDensity" <-
function(node, beg = samplesGetBeg(), end = samplesGetEnd(),
firstChain = samplesGetFirstChain(), lastChain = samplesGetLastChain(),
thin = samplesGetThin(), plot = TRUE, mfrow = c(3, 2), ask = NULL,
ann = TRUE, ...)
# Plot posterior density
{
    if(is.null(ask)) {
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
    if (plot) {
        if (is.R())
            par(mfrow = mfrow, ask = ask, ann = ann)
        else
            par(mfrow = mfrow, ask = ask)
    }
    result <- sapply(mons, plotDensity, plot=plot, ...)
    if (!is.R())
        invisible()
    else {
        if(plot) invisible(result)
        else     return(result)
    }
}
