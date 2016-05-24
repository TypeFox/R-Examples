radialtext <- function(x, center=c(0,0), start=NA, middle=1, end=NA, angle=0, deg=NA,
    expand=0, stretch=1, nice=TRUE, cex=NA, ...) 
    {
    oldcex <- par("cex")
    if (is.na(cex))
        cex <- oldcex
    par(cex=cex)
    if (is.na(deg))
        deg <- angle*180/pi
    deg <- deg %% 360
    chardeg <- deg
    if (nice && (deg > 90 && deg < 270))
        {
        chardeg <- (deg + 180) %% 360
        x <- paste(rev(unlist(strsplit(x, ""))),collapse="")
        }
    angle <- deg*pi/180
    xvec <- strsplit(x, "")[[1]]
    lenx <- length(xvec)
    xwidths <- stretch * strwidth(xvec)
    xwrange <- range(xwidths)
    xwidths[xwidths < xwrange[2]/2] <- xwrange[2]/2 # Make really narrow characters wider
    # Compute expand factor for each character and adjust character widths accordingly.
    chexp <- rep(1, lenx)
    if (!is.na(start) && expand != 0)
        {
        # Note: width of each char changes in succession, affecting sizes AND positions.
        expfactor <- expand/start
        deltar <- 0
        for (xchar in 2:lenx)
            {
            deltar <- deltar + xwidths[xchar-1]
            expansion <- 1+deltar*expfactor
            if (expansion < 0.1) expansion <- 0.1
            chexp[xchar] <- expansion
            xwidths[xchar] <- xwidths[xchar]*expansion
            }
        }
    # Find start distance
    if (is.na(start))
        {
        if (is.na(end))
            start <- middle - sum(xwidths)/2
        else
            start <- end - sum(xwidths)
        }
    cosang <- cos(angle)
    sinang <- sin(angle)
    charstart <- c(start, start + cumsum(xwidths)[-lenx])
    charpos <- charstart + xwidths/2
    xylim <- par("usr")
    plotdim <- par("pin")
    ymult <- (xylim[4] - xylim[3])/(xylim[2] - xylim[1]) * plotdim[1]/plotdim[2]
    for (xchar in 1:lenx)
        {
        par(cex=cex*chexp[xchar])
        text(center[1] + charpos[xchar] * cosang, 
            center[2] + charpos[xchar] * ymult * sinang, xvec[xchar], 
            adj=c(0.5, 0.5), srt=chardeg, ...)
        }
    par(cex=oldcex)
    }
