rp.plot3d <- function (x, y, z, xlab = NA, ylab = NA, zlab = NA, 
    axes = TRUE, new.window = TRUE, type = "p", size = 3, col = "red", 
    xlim = NA, ylim = NA, zlim = NA, plot = TRUE, ...) {
    	
    if (require(rgl)) {
    	
        xname <- deparse(substitute(x))
        if (!missing(y)) yname <- deparse(substitute(y))
        if (!missing(z)) zname <- deparse(substitute(z))
        if (is.data.frame(x)) x <- as.matrix(x)
        
        if (is.matrix(x)) {
           if (!is.null(colnames(x))) xname <- colnames(x)
           if (ncol(x) >= 3) {
              y <- x[ , 2]
              z <- x[ , 3]
              x <- x[ , 1]
           }
           else
              stop("x is a matrix with fewer than three columns.")
           if (is.na(xlab))
              xlab <- if (length(xname) > 1) xname[1] else paste(xname, "1", sep = "-")
           if (is.na(ylab))
              ylab <- if (length(xname) > 1) xname[2] else paste(xname, "2", sep = "-")
           if (is.na(zlab))
              zlab <- if (length(xname) > 1) xname[3] else paste(xname, "3", sep = "-")
        }
        else {
           if (missing(y) | missing(z)) stop("too few arguments.")
           if (is.na(xlab)) xlab <- xname
           if (is.na(ylab)) ylab <- yname
           if (is.na(zlab)) zlab <- zname
        }
        
        xrange <- xlim
        yrange <- ylim
        zrange <- zlim
        ind <- !is.na(x + y + z)
        if (length(col) == length(x))
            ind <- (ind & (!is.na(col)))
        if (!all(ind))
            cat("Warning: missing data removed. \n")
        if (any(is.na(xlim))) {
            xrange[1] <- min(x[ind]) - 0.05 * diff(range(x[ind]))
            xrange[2] <- max(x[ind]) + 0.05 * diff(range(x[ind]))
        }
        if (any(is.na(ylim))) {
            yrange[1] <- min(y[ind]) - 0.05 * diff(range(y[ind]))
            yrange[2] <- max(y[ind]) + 0.05 * diff(range(y[ind]))
        }
        if (any(is.na(zlim))) {
            zrange[1] <- min(z[ind]) - 0.05 * diff(range(z[ind]))
            zrange[2] <- max(z[ind]) + 0.05 * diff(range(z[ind]))
        }
        xscale <- pretty(xrange)
        yscale <- pretty(yrange)
        zscale <- pretty(zrange)
        xscale <- xscale[xscale >= xrange[1] & xscale <= xrange[2]]
        yscale <- yscale[yscale >= yrange[1] & yscale <= yrange[2]]
        zscale <- zscale[zscale >= zrange[1] & zscale <= zrange[2]]
        xadj1 <- mean(xrange)
        yadj1 <- mean(yrange)
        zadj1 <- mean(zrange)
        xadj2 <- diff(xrange)/2
        yadj2 <- diff(yrange)/2
        zadj2 <- diff(zrange)/2
        x.orig <- x
        y.orig <- y
        z.orig <- z
        x <- (x - xadj1)/xadj2
        y <- (y - yadj1)/yadj2
        z <- (z - zadj1)/zadj2
        xscale.adj <- (xscale - xadj1)/xadj2
        yscale.adj <- (yscale - yadj1)/yadj2
        zscale.adj <- (zscale - zadj1)/zadj2
        rx <- c(-1, 1)
        ry <- c(-1, 1)
        rz <- c(-1, 1)
        if (plot) {
           if (new.window) {
              open3d()
              bg3d(col = c("white", "black"))
              }
           else
              rgl.clear()
           rgl.viewpoint(-30, 30, fov = 1)
           if (axes) {
               rgl.lines(rx[c(1, 2, 2, 2, 2, 1, 1, 1)], ry[rep(1,
                   8)], rz[c(1, 1, 1, 2, 2, 2, 2, 1)], col = "black")
               rgl.lines(rx[c(1, 2, 2, 2, 2, 1, 1, 1)], ry[rep(2,
                   8)], rz[c(1, 1, 1, 2, 2, 2, 2, 1)], col = "black")
               for (i in 1:2) for (j in 1:2) rgl.lines(rx[c(i, i)],
                   ry[c(1, 2)], rz[c(j, j)], col = "black")
               rgl.texts(mean(rx), min(rx), min(rx), "")
               delta <- 0.1
               nyticks <- length(yscale)
               if (nyticks/2 - floor(nyticks/2) > 0)
                   ypos <- 1/(nyticks - 1)
               else ypos <- 0
               rgl.texts(c(0, -1 - 2 * delta, -1 - 2 * delta), c(-1 -
                   2 * delta, ypos, -1 - 2 * delta), c(1 + 2 * delta,
                   -1 - 2 * delta, 0), c(xlab, ylab, zlab), adj = c(0.5,
                   0.5), col = "blue")
               rgl.texts((xscale - xadj1)/xadj2, -1 - delta, 1 +
                   delta, as.character(xscale), col = "black")
               rgl.texts(-1 - delta, (yscale - yadj1)/yadj2, -1 -
                   delta, as.character(yscale), col = "black")
               rgl.texts(-1 - delta, -1 - delta, (zscale - zadj1)/zadj2,
                   as.character(zscale), col = "black")
               scaling <- function(x, y, z) list(x = x, y = y, z = z)
               rgl.segments(xscale.adj, -1, 1, xscale.adj, -1 -
                   delta/4, 1 + delta/4, scaling = scaling, col = "black")
               rgl.segments(-1, yscale.adj, -1, -1 - delta/4, yscale.adj,
                   -1 - delta/4, scaling = scaling, col = "black")
               rgl.segments(-1, -1, zscale.adj, -1 - delta/4, -1 -
                   delta/4, zscale.adj, scaling = scaling, col = "black")
               }
           if (!(type == "n")) {
              ind1 <- ((x.orig >= xrange[1]) & (x.orig <= xrange[2]) &
                  (y.orig >= yrange[1]) & (y.orig <= yrange[2]) &
                  (z.orig >= zrange[1]) & (z.orig <= zrange[2]))
              ind <- (ind1 & ind)
              if (length(col) == length(x.orig))
                  clr <- col[ind]
              else clr <- col
              rgl.points(x[ind], y[ind], z[ind], size = size, col = clr)
              }
           }
        scaling <- function(x, y, z) {
            xx <- (x - xadj1)/xadj2
            yy <- (y - yadj1)/yadj2
            zz <- (z - zadj1)/zadj2
            list(x = xx, y = yy, z = zz)
        }
        invisible(scaling)
    }
    else {
        warning("Package rgl is not installed.")
    }
}

rgl.segments <- function(x0, y0, z0, x1, y1, z1, scaling, ...) {
         a <- scaling(c(rbind(x0, x1)), c(rbind(y0, y1)), c(rbind(z0, z1)))
         rgl.lines(a$x, a$y, a$z, ...)
         } 
         
