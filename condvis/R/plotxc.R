plotxc <-
function (xc, xc.cond, name = NULL, select.colour = NULL, select.lwd = NULL, 
    cex.axis = NULL, cex.lab = NULL, tck = NULL, select.cex = 1, ...)
{
    select.colour <- if (is.null(select.colour))
        "black"
    else select.colour
    select.lwd <- if (is.null(select.lwd))
        2
    else select.lwd
    mar <- if (!exists("mar"))
        c(3, 3, 0.5, 0.5)
    else mar
    cex.axis <- if (identical(version$os, "linux-gnu"))
        1
    else if (is.null(cex.axis))
            0.7
        else cex.axis
    cex.lab <- if (identical(version$os, "linux-gnu"))
        1
    else if (is.null(cex.lab))
            0.8
        else cex.lab
    tck <- if (is.null(tck))
        - 0.2
        else tck
    par(mar = mar)
    par(mgp = c(1.5, 0.5, 0))
    if (is.vector(xc) | is.factor(xc)){
        if (!is.factor(xc)){
            histmp <- hist(xc, xlab = name, ylab = "", main = "", cex.axis = 
                cex.axis, cex.lab = cex.lab, tcl = tck, mgp = c(1.5, 0.5, 0.1))
            lines(x = rep(xc.cond, 2L), y = c(0, max(histmp$counts)), col = 
                select.colour, lwd = select.lwd)
            plot.type <- "histogram"
        } else {
            bartmp <- barplot2(table(xc), main = "", xlab = name,
                cex.axis = cex.axis, cex.lab = cex.lab, tcl = tck)
            factorcoords <- data.frame(level = levels(xc),
			    x = - 0.5 + 1.2 * (1:length(levels(xc))))
            barindex <- factorcoords$level == as.character(xc.cond)
            rect(xleft = bartmp$w.l[barindex], xright = bartmp$w.r[barindex],
                ybottom = 0, ytop = bartmp$height[barindex], col = select.colour
                , density = -1)
            plot.type <- "barplot"
            xc.cond <- factor(xc.cond, levels(xc))
        }
    } else {
        if (is.data.frame(xc) & identical(ncol(xc), 2L)){
            are.factors <- vapply(xc,is.factor, logical(1))
            if (all(are.factors)){
                sptmp <- spineplot2(table(xc), ...)
                xmatch <- as.character(xc.cond[, 1]) == levels(xc[, 1])
                ymatch <- as.character(xc.cond[, 2]) == levels(xc[, 2])
                xlev <- levels(xc[, 1])[xmatch]
                ylev <- levels(xc[, 2])[ymatch]
                match.index <- (sptmp$xnames == xlev) & (sptmp$ynames == ylev)
                rect(xleft = sptmp$xleft[match.index],
                    ybottom = sptmp$ybottom[match.index],
                    xright = sptmp$xright[match.index],
                    ytop = sptmp$ytop[match.index],
                    col = select.colour, density = -1)
                axis(1, at = ((sptmp$xat[1L:sptmp$nx] + sptmp$xat[2L:(sptmp$nx
                    + 1L)] - sptmp$off)/2)[xmatch],
                    labels = unique(sptmp$xnames)[xmatch], tick = FALSE,
                    col.axis = select.colour)
                axis(2, at = sptmp$yat[ymatch],
                    labels = unique(sptmp$ynames)[ymatch],
                    col.axis = select.colour, tick = FALSE)
                plot.type <- "spineplot"
            } else {
                if (any(are.factors)){
                    boxx <- xc[, are.factors]
                    boxy <- xc[, !are.factors]
                    boxtmp <- boxplot(boxy ~ boxx, xlab = name[are.factors],
                        ylab = name[!are.factors], cex.axis = cex.axis,
                        cex.lab = cex.lab)
                    factorcoords <- data.frame(
                        level = levels(xc[, are.factors]),
                        x = 1:length(levels(xc[, are.factors])))
                    abline(v = factorcoords$x[as.character(factorcoords$level) 
                        == as.character(xc.cond[,are.factors])], h = xc.cond[
                        !are.factors], lwd = select.lwd, col = select.colour)
                    plot.type <- "boxplot"
                    xc <- xc[, order(!are.factors)]
                    xc.cond <- data.frame(factor(xc.cond[, are.factors], 
                        levels(boxx)), xc.cond[, !are.factors])
                    name <- name[order(!are.factors)]
                    names(xc.cond) <- name
                } else {
                    if (nrow(xc) > 2000 && requireNamespace("gplots", quietly = 
                        TRUE)){
                        b <- seq(0.35, 1, length.out = 16)
                        gplots::hist2d(xc[, 1], xc[, 2], nbins = 50, col = 
                            c("white", rgb(1 - b, 1 - b, 1 - b)), xlab = 
                            colnames(xc)[1], ylab = colnames(xc)[2], cex.axis = 
                            cex.axis, cex.lab = cex.lab, tcl = tck)
                        box()
                    } else {
                        plot.default(xc[, 1], xc[, 2], xlab = colnames(xc)[1],
                            ylab = colnames(xc)[2], cex.axis = cex.axis,
                            cex.lab = cex.lab, tcl = tck, cex = select.cex)
                    }        
                    abline(v = xc.cond[1], h = xc.cond[2], lwd = select.lwd,
                        col = select.colour)
                    plot.type <- "scatterplot"
                }
            }
        } else stop("Unexpected value for 'xc'")
    }
    structure(list(xc = xc, xc.cond.old = xc.cond, name = name, select.colour = 
        select.colour, mar = mar, select.lwd = select.lwd, select.cex = 
        select.cex, cex.axis = cex.axis, cex.lab = cex.lab, tck = tck, device = 
        dev.cur(), usr = par("usr"), screen = screen(), screen.coords = par(
        "fig"), plot.type = plot.type, sptmp = if(exists("sptmp")) sptmp else 
        NULL, factorcoords = if(exists("factorcoords")) factorcoords else NULL, 
        histmp = if(exists("histmp")) histmp else NULL, bartmp = if(exists(
        "bartmp")) bartmp else NULL, boxtmp = if(exists("boxtmp")) boxtmp else 
        NULL, ...), class = "xcplot")
}
