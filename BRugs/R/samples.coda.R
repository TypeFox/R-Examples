"samplesCoda" <- function(node, stem, beg = samplesGetBeg(), 
    end = samplesGetEnd(), firstChain = samplesGetFirstChain(), 
    lastChain = samplesGetLastChain(), thin = samplesGetThin())
{
# Write out CODA files
    if(!is.character(node) || length(node)!=1)
        stop("'node' must be character of length 1")
    if(!is.character(stem) || length(stem)!=1)
        stop("'stem' must be character of length 1")
    if(dirname(stem) == ".") 
        stem <- file.path(getwd(), basename(stem))
    
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
    command <- paste(.SamplesGlobalsCmd(node), ";SamplesEmbed.StatsGuard;",
                     "SamplesEmbed.CODA(", sQuote(stem), ")")
    .CmdInterpreter(command)
    buffer()
}
