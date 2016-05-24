plot.tile.list <- function (x, verbose = FALSE, close = FALSE, pch = 1,
                            fillcol = getCol(x,warn=warn), col.pts=NULL,
                            border=NULL, showpoints = TRUE, add = FALSE,
                            asp = 1, clipp=NULL, xlab = "x", ylab = "y",
                            main = "", warn=FALSE, ...) {
    object <- x
    if (!inherits(object, "tile.list")) 
        stop("Argument \"object\" is not of class tile.list.\n")
    clip  <- !is.null(clipp)
    n     <- length(object)
    rw    <- attr(object, "rw")
    rx    <- rw[1:2]
    ry    <- rw[3:4]
    x.pts <- unlist(lapply(object, function(w) {
        w$pt[1]
    }))
    y.pts <- unlist(lapply(object, function(w) {
        w$pt[2]
    }))
    if (!add) 
        plot(0, 0, type = "n", asp = asp, xlim = rx, ylim = ry, 
            xlab = xlab, ylab = ylab, main = main)
    fillcol <- apply(col2rgb(fillcol, TRUE), 2, function(x) {
        do.call(rgb, as.list(x/255))
    })
    fillcol <- rep(fillcol, length = length(object))
    hexbla <- do.call(rgb, as.list(col2rgb("black", TRUE)/255))
    hexwhi <- do.call(rgb, as.list(col2rgb("white", TRUE)/255))
    if(is.null(col.pts)){
        col.pts <- ifelse(fillcol == hexbla, hexwhi, hexbla)
    } else {
        col.pts <- apply(col2rgb(col.pts, TRUE), 2, function(x) {
            do.call(rgb, as.list(x/255))
        })
        col.pts <- rep(col.pts, length = length(object))
    }
    if(is.null(border))
        border <- if(all(fillcol == hexbla)) hexwhi else hexbla
    else if(length(border) > 1)
        stop("Argument \"border\" must be a scalar or NULL.\n")
    lnwid <- if(all(fillcol == hexbla)) 2 else 1
    pch <- rep(pch,n)
    okn <- logical(n)
    for(i in 1:n) {
        if(clip) {
            if(requireNamespace("polyclip",quietly=TRUE)) {
                pgon <- polyclip::polyclip(object[[i]],clipp)
                ok   <- length(pgon) > 0
            } else {
                stop("Cannot clip the tiles; package \"polyclip\" not available.\n")
            }
        } else {
            pgon <- list(object[[i]])
            ok <- TRUE
        }
        okn[i] <- ok
        inner <- !any(object[[i]]$bp)
        for(ii in seq(along=pgon)){
            ptmp <- pgon[[ii]]
            polygon(ptmp,col=fillcol[i],border=NA)
            if (close | inner) { 
                polygon(ptmp,col = NA, border = border, lwd = lnwid)
            } else {
                x <- ptmp$x
                y <- ptmp$y
                ni <- length(x)
                for (j in 1:ni) {
                    jnext <- if (j < ni) j + 1 else 1
                    do.it <- mid.in(x[c(j, jnext)], y[c(j, jnext)], rx, ry)
                    if (do.it) 
                        segments(x[j], y[j], x[jnext], y[jnext],
                                 col = border, lwd = lnwid)
                }
            }
         }
         if(ok & verbose) {
             if(showpoints) points(object[[i]]$pt[1], object[[i]]$pt[2],
                                   pch = pch[i], col = col.pts[i])
             if(i < n) readline(paste("i = ",i,"; Go? ",sep=""))
             if(i == n) cat("i = ",i,"\n",sep="")
         }
    }
    if (showpoints & !verbose) 
    points(x.pts[okn], y.pts[okn], pch = pch[okn], col = col.pts[okn])
    invisible()
}
