drawmap <- function(data, map,
                    regionvar=2, plotvar=3,
                    limits,
                    cols="hcl", nrcolors=100, swapcolors=FALSE,
                    pcat=FALSE,
                    hcl.par=list(h=c(120, 0), c=60, l=c(45,60), power=1.2),
                    hsv.par=list(s=1, v=1), legend=TRUE, drawnames=FALSE, cex.names=0.7, 
                    cex.legend=0.7, mar.min=2, density=15, ...)
{
    if(! inherits(map,"bnd"))
        stop("Argument 'map' is not an object of class 'bnd'!")

    regions <- names(map)
    S <- length(regions)
    is.in <- attr(map, "is.in")
    height2width <- attr(map, "height2width")
    height2width <- height2width*1.1

    ## draw inner regions last: reorder therefore the map
    surrounding <- attr(map, "surrounding")
    innerRegionInds <- which(sapply(surrounding, length) > 0L)
    if(length(innerRegionInds)){
        regions <- c(regions[- innerRegionInds], regions[innerRegionInds])
        map <- c(map[- innerRegionInds], map[innerRegionInds])
    }

    if(!is.null(mar.min)){
        if(height2width > 1){
            side <- 17.5*(1-1/height2width)+mar.min/height2width
            par(mar=c(mar.min, side, mar.min, side))
        }
        else{
            top <- 17.5*(1-height2width)+mar.min*height2width
            par(mar=c(top,mar.min,top,mar.min))
        }
    }

    black <- grey(0)
    white <- grey(1)

    xmin <- 1:S
    xmax <- 1:S
    ymin <- 1:S
    ymax <- 1:S
    for (i in 1:S) {
        xmin[i] <- min(map[[i]][, 1], na.rm = TRUE)
        xmax[i] <- max(map[[i]][, 1], na.rm = TRUE)
        ymin[i] <- min(map[[i]][, 2], na.rm = TRUE)
        ymax[i] <- max(map[[i]][, 2], na.rm = TRUE)
    }
    xlimits <- c(min(xmin), max(xmax))
    ylimits <- c(min(ymin) - (max(ymax) - min(ymin)) * 0.1, max(ymax))

    if(missing(data)){
        plot(xlimits, ylimits, type = "n", axes = FALSE, xlab="", ylab="", ...)
        for (k in 1:S)
            polygon(map[[k]][, 1], map[[k]][,2], lwd = 0.3, border = black)
    }

    else{
        if(!is.data.frame(data))
            data <- read.table(data, header = TRUE)

        ord <- order(data[, regionvar])
        plotvar <- data[, plotvar]
        plotvar <- plotvar[ord]
        regionvar <- data[, regionvar]
        regionvar <- regionvar[ord]

        if(cols != "hcl" && cols != "hsv" && cols != "grey") {
            nrcolors <- length(cols)
            if (swapcolors == TRUE)
                cols <- rev(cols)
        }
        else{
            if(cols == "hcl"){
                h <- hcl.par$h
                c <- hcl.par$c
                l <- hcl.par$l
                power <- hcl.par$power
            }
            if(cols == "hsv"){
                s <- hsv.par$s
                v <- hsv.par$v
            }
        }


        maxim <- max(plotvar, na.rm = TRUE)
        minim <- min(plotvar, na.rm = TRUE)

        if(cols != "hcl" && cols != "hsv" && cols != "grey") {
            upperlimit <- 1
            lowerlimit <- -1
        }
    	else {
            if (missing(limits)) {
                lowerlimit <- minim
                upperlimit <- maxim
            }
            else {
                lowerlimit <- limits[1]
                upperlimit <- limits[2]

                if (lowerlimit > minim) {
                    plotvar[plotvar < lowerlimit] <- lowerlimit
                    cat(paste("Note: lowerlimit is above minimum value (", lowerlimit, " > ", minim, ")\n"))
                }
                if (upperlimit < maxim) {
                    plotvar[plotvar > upperlimit] <- upperlimit
                    cat(paste("Note: upperlimit is below maximum value (", upperlimit, " < ", maxim, ")\n"))
                }
            }
        }

        if(pcat) {
            nrcolors <- 3
            upperlimit <- 1
            lowerlimit <- -1
            if(cols != "hcl" && cols != "hsv" && cols != "grey")
                cols <- c(cols[1], cols[round(length(cols)/2 + 0.5)], cols[length(cols)])
        }


        fill.colors <- cut(c(lowerlimit, plotvar, upperlimit),nrcolors)
        fill.colors <- fill.colors[c(-1, -length(fill.colors))]
        fill.colors <- as.vector(fill.colors, mode = "numeric")

     	if(cols != "hcl" && cols != "hsv" && cols != "grey"){
            fill.colors <- cols[fill.colors]
            legend.colors <- cols
        }
        else{
            if (cols == "hcl") {
                if (swapcolors == TRUE)
                    h <- rev(h)
                fill.colors <- colorspace::diverge_hcl(nrcolors,h=h,c=c,l=l,power=power)[fill.colors]
                legend.colors <- colorspace::diverge_hcl(nrcolors,h=h,c=c,l=l,power=power)
            }
            if (cols == "hsv") {
                fill.colors <- (fill.colors-1)/(3*(nrcolors-1))
                if (swapcolors == FALSE)
                    fill.colors <- 1/3 - fill.colors
                fill.colors <- hsv(h = fill.colors, s=s, v=v)
                legend.colors <- hsv(h = (0:(nrcolors-1))/(3*(nrcolors-1)), s=s, v=v)
                if (swapcolors == FALSE)
                    legend.colors <- rev(legend.colors)
            }
            if (cols == "grey") {
                fill.colors <- (fill.colors-1)/(nrcolors-1)
                if (swapcolors == TRUE)
                    fill.colors <- 1 - fill.colors
                fill.colors <- grey(fill.colors)
                legend.colors <- grey((0:(nrcolors-1))/(nrcolors-1))
                if (swapcolors == TRUE)
                    legend.colors <- rev(legend.colors)
            }
        }


        plot(xlimits, ylimits, type = "n", axes = FALSE, col = white, xlab="", ylab="", ...)

        if(sum(!is.na(match(regions, regionvar))) == 0)
            warning("map probably doesn't match datafile")
        block1 <- c()
        block2 <- c()
        for (k in 1:S) {
            if (is.na(map[[k]][1, 1]) && is.na(map[[k]][1, 2]))
                block2 <- c(block2, k)
            else
                block1 <- c(block1, k)
        }
        m <- match(regions, regionvar)
        for (k in block1) {
            if (is.na(m[k])) {
                polygon(map[[k]][, 1], map[[k]][, 2], col = white, border = FALSE)
                polygon(map[[k]][, 1], map[[k]][, 2], density  = density, lwd = 0.3, col = black)
            }
            else
                polygon(map[[k]][, 1], map[[k]][, 2], col = fill.colors[m[k]], border = black)
        }
        for (k in block2) {
            if (is.na(m[k])) {
                polygon(map[[k]][-1, 1], map[[k]][-1, 2], col = white, border = FALSE)
                polygon(map[[k]][-1, 1], map[[k]][-1, 2], density = density, lwd = 0.3, col = black)
            }
            else
                polygon(map[[k]][-1, 1], map[[k]][-1, 2], col = fill.colors[m[k]], border = black)
        }

        if (legend == TRUE) {
            ylo <- yro <- ylimits[1] + (0.7 * (ylimits[2] - ylimits[1]))/11
            ylu <- yru <- ylimits[1] + (0.3 * (ylimits[2] - ylimits[1]))/11
            tylu <- tyru <- ylimits[1]
            xlu <- xlo <- xlimits[1] + 0.1 * (xlimits[2] - xlimits[1])
            xru <- xro <- xlimits[1] + 0.4 * (xlimits[2] - xlimits[1])
            step <- (xru - xlu)/nrcolors
            for (i in 0:(nrcolors - 1)) {
                polygon(c(xlo + step * i, xlo + step * (i + 1), xlu + step * (i + 1), xlu + step * i),
                        c(ylo, yro, yru, ylu), col = legend.colors[i + 1], border = legend.colors[i + 1])
            }
            lines(c(xlo, xro, xru, xlu, xlo), c(ylo, yro, yru, ylu, ylo), col = black)
            text(xlu + 0.5 * step, tylu, round(lowerlimit,4), col = black, cex = cex.legend)
            text(xru - 0.5 * step, tyru, round(upperlimit,4), col = black, cex = cex.legend)
            
            if (lowerlimit + (upperlimit - lowerlimit)/3 < 0 && 0 < upperlimit - (upperlimit - lowerlimit)/3) {
                help <- cut(c(0, lowerlimit, upperlimit), nrcolors)
                help <- as.vector(help, mode = "numeric")
                if(nrcolors%%2 == 0)
                    text(xlu + step * (help[1]), tylu, "0", col = black, cex = cex.legend)              
                else
                    text(xlu + step * (help[1] - 0.5), tylu, "0", col = black, cex = cex.legend)
            }
        }
    }
    
    if (drawnames == TRUE) {
        xpos <- (xmin + xmax)/2
        ypos <- (ymin + ymax)/2
        text(xpos, ypos, labels = regions, col = black, cex = cex.names)
    }

    return(invisible())
}

