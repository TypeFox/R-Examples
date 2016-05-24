polys.map <- function (lm, z, scheme = "gray", lab = "", zlim, rev.col = TRUE, ...){

poly2 <- getFromNamespace("poly2", "mgcv") 

    for (i in 1:length(lm)) {
        yr <- range(lm[[i]][, 2], na.rm = TRUE)
        xr <- range(lm[[i]][, 1], na.rm = TRUE)
        if (i == 1) {
            ylim <- yr
            xlim <- xr
        }
        else {
            if (yr[1] < ylim[1]) 
                ylim[1] <- yr[1]
            if (yr[2] > ylim[2]) 
                ylim[2] <- yr[2]
            if (xr[1] < xlim[1]) 
                xlim[1] <- xr[1]
            if (xr[2] > xlim[2]) 
                xlim[2] <- xr[2]
        }
    }
    mar <- par("mar")
    oldpar <- par(mar = c(2, mar[2], 2, 1))

        nz <- names(z)
        nlm <- names(lm)
        if (!is.null(nz) && !is.null(nlm)) {
            if (all.equal(sort(nz), sort(nlm)) != TRUE) 
                stop("names of z and lm must match")
            z <- z[nlm]
        }
        xmin <- xlim[1]
        xlim[1] <- xlim[1] - 0.1 * (xlim[2] - xlim[1])
        n.col <- 100
        if(scheme == "heat")    schem <- heat.colors(n.col + 1)
        if(scheme == "terrain") schem <- terrain.colors(n.col + 1)
        if(scheme == "topo")    schem <- topo.colors(n.col + 1)
        if(scheme == "cm")      schem <- cm.colors(n.col + 1)
        if(scheme == "gray")    schem <- gray(0:n.col/n.col)
        
        if(rev.col == TRUE)     schem <- rev(schem)
        
        if(missing(zlim)) zlim <- range(pretty(z))
        
        for (i in 1:length(lm)) lm[[i]][, 2] <- zlim[1] + (zlim[2] - 
            zlim[1]) * (lm[[i]][, 2] - ylim[1])/(ylim[2] - ylim[1])
        ylim <- zlim
        plot(0, 0, ylim = ylim, xlim = xlim, type = "n", xaxt = "n", 
            bty = "n", xlab = "", ylab = lab, ...)
        for (i in 1:length(lm)) {
            coli <- round((z[i] - zlim[1])/(zlim[2] - zlim[1]) * 
                n.col) + 1
            poly2(lm[[i]], col = schem[coli])
        
        xmin <- min(c(axTicks(1), xlim[1]))
        dx <- (xlim[2] - xlim[1]) * 0.05
        x0 <- xmin - 2 * dx
        x1 <- xmin + dx
        dy <- (ylim[2] - ylim[1])/n.col
        poly <- matrix(c(x0, x0, x1, x1, ylim[1], ylim[1] + dy, 
            ylim[1] + dy, ylim[1]), 4, 2)
        for (i in 1:n.col) {
            polygon(poly, col = schem[i], border = NA)
            poly[, 2] <- poly[, 2] + dy
        }
        poly <- matrix(c(x0, x0, x1, x1, ylim[1], ylim[2], ylim[2], 
            ylim[1]), 4, 2)
        polygon(poly, border = "black")
    }
    par(oldpar)
}




























