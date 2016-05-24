### Plot bgr diagnostic for single component of OpenBUGS name
"plotBgr" <-
    function(node, plot = TRUE, main = NULL, xlab = "iteration", ylab = "bgr", 
             col = c("red", "blue", "green"), bins = 50, ...)
{
    sM <- samplesMonitors(node)
    if(length(sM) > 1 || sM != node)
        stop("node must be a scalar variable from the model, for arrays use samplesBgr")
    if (any(grep("^inference can not be made", sM))) { stop(sM) }
    grid <- bgrGrid(node, bins = bins)

    ## Use a single API call instead of looping API calls over
    ## iterations - more efficient with the Linux helper.
    
    ## find size of available sample at each grid point 
    res <- .OpenBUGS(cmds = c(.SamplesGlobalsCmd(node),
                     as.vector(rbind(paste("SamplesEmbed.end := ", grid, ";"), "SamplesEmbed.SampleSize;"))),
                     cmdtypes = c("CmdInterpreter", rep(c("CmdInterpreter","Integer"), bins)),
                     args=as.list(c(NA, rep(c(NA, NA), bins)))
                     )
    
    args <- list(NA)
    for (i in seq(length=bins)){
        args[[2*i]] <- NA
        args[[2*i + 1]] <- double(res[[2*i + 1]])
    }

    ## get available sample at each grid point 
    res <- .OpenBUGS(cmds =
                     c(.SamplesGlobalsCmd(node),
                       as.vector(rbind(paste("SamplesEmbed.end := ", grid, ";"), "SamplesEmbed.SampleValues;"))), 
                     cmdtypes = c("CmdInterpreter", rep(c("CmdInterpreter","RealArray"), bins)),
                     args=args)

    ## remove junk elements of list, leaving a list of samples for each grid point
    res[c(1, 2*seq(length=bins))] <- NULL

    ## calculate between, within and ratio statistics for each grid point
    bgr <- rbind(grid, sapply(res, bgrPoint))

    yRange <- range(bgr[4,])
    yRange <- c(0, max(c(1.2, yRange[2])))
    nRange <- range(bgr[2,])
    nRange <- c(min(c(0, nRange[1])), nRange[2])
    nDelta <- nRange[2] - nRange[1]
    dRange <- range(bgr[3,])
    dRange <- c(min(c(0, dRange[1])), dRange[2])
    dDelta <- dRange[2] - dRange[1]
    max <- 2 * max(c(nDelta, dDelta))
    bgr[2,] <- bgr[2,] / max
    bgr[3,] <- bgr[3,] / max
    if(plot){
        plot(grid, bgr[4,], ylim = yRange, type = "l", 
             main = if(is.null(main)) node else main, xlab = xlab, ylab = ylab, col = col[1], ...)
        lines(grid, bgr[2,], col = col[2], ...)
        lines(grid, bgr[3,], col = col[3], ...)
    }
    bgr <- data.frame(t(bgr))
    names(bgr) <- c("Iteration", "pooledChain80pct", "withinChain80pct", "bgrRatio")
    bgr$Iteration <- as.integer(bgr$Iteration)
    if(plot) invisible(bgr)
    else return(bgr)

}
