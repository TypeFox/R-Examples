hilbert <- function(xt)
{     
    if(is.ts(xt))
        xt <- as.numeric(xt) 
        
    ndata <- length(xt)
    h <- rep(0, ndata)
    if(ndata %% 2 == 0) {
        h[c(1, ndata/2+1)] <- 1
        h[2:(ndata/2)] <- 2
    } else { 
        h[1] <- 1
        h[2:((ndata+1)/2)] <- 2
    }
    xt <- fft(h * fft(xt), inverse = TRUE) / ndata
    return(xt)
}

unwrap <- function(phase, tol = pi)
{

    d = c(0, -diff(phase))
    p = 2 * pi  * ((d > tol) - (d < -tol))

    unphase = phase + cumsum (p)
    return(unphase)
}

hilbertspec <- function (xt, tt=NULL, central=FALSE)
{
    if (is.null(tt))
        tt <- 1:nrow(xt)
    deltaT <- diff(tt)
    ndata <- nrow(xt)
    nts <- ncol(xt)
    amplitude <- instantfreq <- energy <- NULL
    for (i in 1:nts) {
        tmpxt <- xt[, i]
        htmpxt <- hilbert(tmpxt)
        amplitude <- cbind(amplitude, abs(htmpxt))
        if (central) { #DANIEL BOWMAN and KEEHOON KIM'S EDITS, MARCH 1 2013
            yt <- Im(htmpxt)
            dyt <- c((yt[2]-yt[1])/(tt[2]-tt[1]), diff(yt, 2)/diff(tt, 2), (yt[ndata]-yt[ndata-1])/(tt[ndata]-tt[ndata-1]))
            dxt <- c((tmpxt[2]-tmpxt[1])/(tt[2]-tt[1]), diff(tmpxt, 2)/diff(tt, 2), (tmpxt[ndata]-tmpxt[ndata-1])/(tt[ndata]-tt[ndata-1]))
            tmpfreq <- (xt*dyt - yt*dxt)/((xt^2) + (yt^2))
        } else { #DANIEL BOWMAN and KEEHOON KIM'S EDITS, MARCH 1 2013
            phase <- atan2(Im(htmpxt), Re(htmpxt))
            unphase <- unwrap(phase)
            tmpfreq <- diff(unphase)/deltaT
            tmpfreq <- abs(tmpfreq[-length(tmpfreq)] + tmpfreq[-1])/2
            tmpfreq <- c(tmpfreq[1], tmpfreq, tmpfreq[length(tmpfreq)])
        } 
        instantfreq <- cbind(instantfreq, tmpfreq)
        energy <- c(energy, sum(apply(as.matrix(amplitude^2), 1, sum)))
    }
    energy <- energy/sum(apply(amplitude^2, 1, sum))
    list(amplitude = amplitude, instantfreq = instantfreq/(2 * pi), energy = energy)
}

    
spectrogram <- function(amplitude, freq, tt=NULL, multi=FALSE, nlevel=NULL, size=NULL) {

    amplitude <- as.matrix(amplitude)
    freq <- as.matrix(freq)
    if(!all(dim(amplitude) == dim(freq))) stop("Dimension of amplitude and frequency must be the same.")
     
    ntime <- nrow(amplitude)
    nseries <- ncol(amplitude)
    
    if(is.null(tt)) tt <- 1:ntime
    
    nlevel <- min(nlevel, nrow(freq))
    nnrow <- min(size[1], ntime)
    nncol <- min(size[2], nrow(freq))
    
    totamp <- totfreq <- NULL
    for (i in 1:nseries) {
        totamp <- rbind(totamp, amplitude[, i])
        totfreq <- rbind(totfreq, cbind(tt, freq[, i]))
    }   
    rangeamp <- range(totamp)
    rangefreq <- range(totfreq[, 2])
         
    if(!multi) {
        tmpim <- as.image(c(totamp), x = totfreq, nrow = nnrow, ncol = nncol)
        #tmpim$z[is.na(tmpim$z)] <- 0
        #tmpim$z <- scale(tmpim$z)
        image.plot.ts(tmpim, tt = tt, col = gray(nlevel:0/nlevel), ylim = rangefreq, zlim = rangeamp)  
        box()
    } else {
        if (nseries != 1) par(mfrow=c(nseries %/% 2 + nseries %% 2, 2), mar = c(4, 4, 1, 4) + 0.1)
        for (i in 1:nseries) {
            tmpim <- as.image(amplitude[, i], x = cbind(tt, freq[, i]), nrow = nnrow, ncol = nncol)
            #tmpim$z <- scale(tmpim$z)
            image.plot.ts(tmpim, tt=tt, col = gray(nlevel:0/nlevel), ylim = rangefreq, zlim = range(amplitude[, i]))
            box()
        }    
    }
}  



image.plot.ts <- function (..., tt, add = FALSE, nlevel = 64, legend.shrink = 0.9, 
    legend.width = 1.2, legend.mar = NULL, graphics.reset = FALSE, 
    horizontal = FALSE, bigplot = NULL, smallplot = NULL, legend.only = FALSE, 
    col = tim.colors(nlevel)) 
{
    old.par <- par(no.readonly = TRUE)
    info <- imageplot.info(...)
    if (add) {
        big.plot <- old.par$plt
    }
    if (legend.only) {
        graphics.reset <- TRUE
    }
    if (is.null(legend.mar)) {
        legend.mar <- ifelse(horizontal, 3.1, 5.1)
    }
    temp <- imageplot.setup(add = add, legend.shrink = legend.shrink, 
        legend.width = legend.width, legend.mar = legend.mar, 
        horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
    smallplot <- temp$smallplot
    bigplot <- temp$bigplot
    if (!legend.only) {
        if (!add) {
            par(plt = bigplot)
        }
        if (inherits(tt, "POSIXt")) {
            image(..., xaxt = "n", add = add, col = col)
            axis.POSIXct(1, tt)
        }
        else
            image(..., add = add, col = col)       
        big.par <- par(no.readonly = TRUE)
    }
    if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
        par(old.par)
        stop("plot region too small to add legend\n")
    }
    ix <- 1
    minz <- info$zlim[1]
    maxz <- info$zlim[2]
    binwidth <- (maxz - minz)/nlevel
    midpoints <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
    iy <- midpoints
    iz <- matrix(iy, nrow = 1, ncol = length(iy))
    breaks <- list(...)$breaks
    if (!horizontal) {
        par(new = TRUE, pty = "m", plt = smallplot, err = -1)
        if (is.null(breaks)) {
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
                ylab = "", col = col)
        }
        else {
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
                ylab = "", col = col, breaks = breaks)
        }
        axis(4, mgp = c(3, 1, 0), las = 2)
        box()
    }
    else {
        par(new = TRUE, pty = "m", plt = smallplot, err = -1)
        if (is.null(breaks)) {
            image(iy, ix, t(iz), yaxt = "n", xlab = "", ylab = "", 
                col = col)
        }
        else {
            image(iy, ix, t(iz), yaxt = "n", xlab = "", ylab = "", 
                col = col, breaks = breaks)
        }
        box()
    }
    mfg.save <- par()$mfg
    if (graphics.reset | add) {
        par(old.par)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
    else {
        par(big.par)
        par(plt = big.par$plt, xpd = FALSE)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
}
