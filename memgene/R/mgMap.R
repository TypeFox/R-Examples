mgMap <-
function(coords, memgene, wid=NULL, hei=NULL, dev.open=FALSE, add.plot=FALSE, legend=FALSE, ...) {
    
    ## FUNCTIONS FROM NUMERICAL ECOLOGY WITH R (Borcard et al, 2011)
    ## included here to support MEM analyses
    sr.value <- function (dfxy, z, xax = 1, yax = 2, method = c("bubble",
        "greylevel"), zmax = NULL, csize = 1, cpoint = 0, pch = 20,
        clegend = 0.75, neig = NULL, cneig = 1, xlim = NULL, ylim = NULL,
        grid = TRUE, addaxes = TRUE, cgrid = 0.75, include.origin = TRUE,
        origin = c(0, 0), sub = "", csub = 1, possub = "topleft",
        pixmap = NULL, contour = NULL, area = NULL, add.plot = FALSE) {
    
    #
    # Slightly modified version of ade4's s.value() graphical function.
    # Draws round instead of square bubbles in some plots when argument 
    # "bubble" is called.
    #
    # License: GPL-2
    # Author of the original function s.value: Daniel Chessel
    # Modification: Francois Gillet, 25 August 2012
    #
        dfxy <- data.frame(dfxy)
        if (length(z) != nrow(dfxy))
            stop(paste("Non equal row numbers", nrow(dfxy), length(z)))
        opar <- par(mar = par("mar"))
        on.exit(par(opar))
        par(mar = c(0.1, 0.1, 0.1, 0.1))
        coo <- scatterutil.base(dfxy = dfxy, xax = xax, yax = yax,
            xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes,
            cgrid = cgrid, include.origin = include.origin, origin = origin,
            sub = sub, csub = csub, possub = possub, pixmap = pixmap,
            contour = contour, area = area, add.plot = add.plot)
        if (!is.null(neig))
        {
            if (is.null(class(neig))) neig <- NULL
            if (class(neig) != "neig") neig <- NULL
            deg <- attr(neig, "degrees")
            if (length(deg) != length(coo$x)) neig <- NULL
        }
        if (!is.null(neig))
        {
            fun <- function(x, coo)
            {
                segments(coo$x[x[1]], coo$y[x[1]], coo$x[x[2]], coo$y[x[2]],
                    lwd = par("lwd") * cneig)
            }
            apply(unclass(neig), 1, fun, coo = coo)
        }
        method <- method[1]
        if (method == "greylevel")
        {
            br0 <- pretty(z, 6)
            nborn <- length(br0)
            coeff <- diff(par("usr")[1:2])/15
            numclass <- cut.default(z, br0, include.lowest = TRUE, labels = FALSE)
            valgris <- seq(1, 0, le = (nborn - 1))
            h <- csize * coeff
            for (i in 1:(nrow(dfxy)))
                {
                    symbols(coo$x[i], coo$y[i], circles = h/2, 
                        bg = gray(valgris[numclass[i]]),
                        add = TRUE, inches = FALSE)
                }
            scatterutil.legend.circle.grey(br0, valgris, h/2, clegend)
            if (cpoint > 0) points(coo$x, coo$y, pch = pch, cex = par("cex") * cpoint)
        }
        else if (method == "bubble")
        {
            coeff <- diff(par("usr")[1:2])/15
            sq <- sqrt(abs(z))
            if (is.null(zmax)) zmax <- max(abs(z))
            w1 <- sqrt(zmax)
            sq <- csize * coeff * sq/w1
            for (i in 1:(nrow(dfxy)))
            {
                if (sign(z[i]) >= 0)
                {
                    symbols(coo$x[i], coo$y[i], circles = sq[i]/2, bg = "black", 
                        fg = "white", add = TRUE, inches = FALSE)
                }
                else
                {
                    symbols(coo$x[i], coo$y[i], circles = sq[i]/2, bg = "white", 
                        fg = "black", add = TRUE, inches = FALSE)
                }
            }
            br0 <- pretty(z, 4)
            l0 <- length(br0)
            br0 <- (br0[1:(l0 - 1)] + br0[2:l0])/2
            sq0 <- sqrt(abs(br0))
            sq0 <- csize * coeff * sq0/w1
            sig0 <- sign(br0)
            if (clegend > 0) scatterutil.legend.bw.circle(br0, sq0, sig0, clegend)
            if (cpoint > 0) points(coo$x, coo$y, pch = pch, cex = par("cex") * cpoint)
        }
        else if (method == "circlesize") print("not yet implemented")
        if (!add.plot) box()
        invisible(match.call())
    }
    
    
    
    scatterutil.legend.bw.circle <- function (br0, sq0, sig0, clegend) {
        br0 <- round(br0, digits = 6)
        cha <- as.character(br0[1])
        for (i in (2:(length(br0)))) cha <- paste(cha, br0[i], sep = " ")
        cex0 <- par("cex") * clegend
        yh <- max(c(strheight(cha, cex = cex0), sq0))
        h <- strheight(cha, cex = cex0)
        y0 <- par("usr")[3] + yh/2 + h/2
        ltot <- strwidth(cha, cex = cex0) + sum(sq0) + h
        rect(par("usr")[1] + h/4, y0 - yh/2 - h/4, 
            par("usr")[1] + ltot + h/4, y0 + yh/2 + h/4, col = "white")
        x0 <- par("usr")[1] + h/2
        for (i in (1:(length(sq0))))
        {
            cha <- br0[i]
            cha <- paste(" ", cha, sep = "")
            xh <- strwidth(cha, cex = cex0)
            text(x0 + xh/2, y0, cha, cex = cex0)
            z0 <- sq0[i]
            x0 <- x0 + xh + z0/2
            if (sig0[i] >= 0)
                symbols(x0, y0, circles = z0/2, bg = "black", fg = "white",
                    add = TRUE, inches = FALSE)
            else symbols(x0, y0, circles = z0/2, bg = "white", fg = "black",
                add = TRUE, inches = FALSE)
            x0 <- x0 + z0/2
        }
        invisible()
    }
    
    
    
    scatterutil.legend.circle.grey <- function (br0, valgris, h, clegend) {
        if (clegend <= 0) return(invisible())
        br0 <- round(br0, digits = 6)
        nborn <- length(br0)
        cex0 <- par("cex") * clegend
        x0 <- par("usr")[1] + h
        x1 <- x0
        for (i in (2:(nborn)))
        {
            x1 <- x1 + h
            cha <- br0[i]
            cha <- paste(cha, "]", sep = "")
            xh <- strwidth(cha, cex = cex0)
            if (i == (nborn)) break
            x1 <- x1 + xh + h
        }
        yh <- max(strheight(paste(br0), cex = cex0), h)
        y0 <- par("usr")[3] + yh/2 + h/2
        rect(par("usr")[1] + h/4, y0 - yh/2 - h/4, x1 - h/4, y0 + yh/2 + h/4, 
            col = "white")
        x0 <- par("usr")[1] + h
        for (i in (2:(nborn)))
        {
            symbols(x0, y0, circles = h/2, bg = gray(valgris[i - 1]), add = TRUE, 
                inches = FALSE)
            x0 <- x0 + h
            cha <- br0[i]
            if (cha < 1e-05) cha <- round(cha, digits = 3)
            cha <- paste(cha, "]", sep = "")
            xh <- strwidth(cha, cex = cex0)
            if (i == (nborn)) break
            text(x0 + xh/2, y0, cha, cex = cex0)
            x0 <- x0 + xh + h
        }
        invisible()
    }
    
    
    
    
    
    ## mgMap specific code
    
    memgene <- as.matrix(memgene)
    
    if (!legend) {
        clegend <- -1
    }
    else {
        clegend <- 0.75
    }
    
    if (is.null(wid) && is.null(hei) && (ncol(memgene)==2) && !add.plot) {
        sideBySide <- TRUE
    }
    else {
        sideBySide <- FALSE
    }

    if (sideBySide) {
        if (!dev.open) {
            dev.new(width=9, height=4.5)
        }
        par(mfcol=c(1,2))
        par(mar=c(2,2,2,2))
        plot(coords, type="n", main="", xlab="", ylab="")
        sr.value(coords, memgene[,1], add.plot=TRUE, clegend=clegend, ...)
        plot(coords, type="n", main="", xlab="", ylab="")
        sr.value(coords, memgene[,2], add.plot=TRUE, clegend=clegend, ...)
    }
    else {
        if (is.null(wid)) {
            wid <- 7
            hei <- 7
        }
        for (i in 1:ncol(memgene)) {
            if (!add.plot) {
                if (!dev.open) {
                    dev.new(width=wid, height=hei)
                }
                par(mar=c(2,2,2,2))
                plot(coords, type="n", main="", xlab="", ylab="")
                sr.value(coords, memgene[,i], add.plot=TRUE, clegend=clegend, ...)

            }
            else {
                sr.value(coords, memgene[,i], add.plot=TRUE, clegend=clegend, ...)
            }
        }
    }
}
