`plot.predcoca` <- function(x,
                            which = c("response", "predictor"),
                            choices = 1:2,
                            display = c("species", "sites"),
                            type,
                            xlim = NULL,
                            ylim = NULL,
                            main = "", sub = "",
                            ylab, xlab,
                            ann = par("ann"),
                            axes = TRUE,
                            ...) {
    ## process the scores to display
    if(missing(display))
        display <-  c("species", "sites")
    display <- match.arg(display, several.ok = TRUE)
    ## what are we plotting, response or predictor?
    which <- match.arg(which)
    ## and map to X and Y for extraction
    WHICH <- ifelse(which == "response", "Y", "X")
    ## need two and only two axes to plot
    if(length(choices) != 2)
        stop("Exactly two axes should be specified in `choices`")
    ## extract the scores
    scrs <- scores(x, choices = choices, display = display)
    ## then extract the response or predictor scores
    scrs <- lapply(scrs, `[[`, WHICH)
    ## what type of plot?
    TYPES <- c("text", "points", "none")
    if (missing(type)) { ## work out whether to plot by text or points
        nitlimit <- 80
        nit <- max(nrow(scrs$species), nrow(scrs$sites))
        type <- if (nit > nitlimit)
            "points"
        else "text"
    } else type <- match.arg(type, TYPES)
    ## compute xy coords for each set of scores
    xy <- lapply(scrs, xy.coords)
    ## process axis limits if non supplied
    if (is.null(xlim))
        xlim <- range(sapply(xy, function(x) range(x$x[is.finite(x$x)])))
    if (is.null(ylim))
        ylim <- range(sapply(xy, function(x) range(x$y[is.finite(x$y)])))
    ## process x/y labels
    if(missing(xlab)) {
        xlabs <- sapply(xy, `[[`, "xlab")
        xlab <- xlabs[!is.null(xlabs)][1]
        if(!is.null(x$lambda)) {
            eigx <- round(x$lambda[choices[1]], 4)
            xlab <- bquote(.(xlab) ~~ (lambda[.(choices[1])] == .(eigx)))
        } else {
            xlab <- bquote(.(xlab))
        }
    }
    if(missing(ylab)) {
        ylabs <- sapply(xy, `[[`, "ylab")
        ylab <- ylabs[!is.null(ylabs)][1]
        if(!is.null(x$lambda)) {
            eigy <- round(x$lambda[choices[2]], 4)
            ylab <- bquote(.(ylab) ~~ (lambda[.(choices[2])] == .(eigy)))
        } else {
            ylab <- bquote(.(ylab))
        }
    }
    #opar <- par(no.readonly=TRUE)
    #on.exit(par(opar))
    ## plotting
    plot.new()
    plot.window(xlim, ylim, ...)
    abline(h = 0, lty = "dashed", col = "grey")
    abline(v = 0, lty = "dashed", col = "grey")
    if(!is.null(scrs$species)) {
        if(type == "text") {
            text(scrs$species, rownames(scrs$species), col = "red",
                 cex = 0.7, ...)
        }
        if(type == "points"){
            points(scrs$species, col = "red", pch = 3, cex = 0.7, ...)
        }
    }
    if(!is.null(scrs$sites)) {
        if(type == "text") {
            text(scrs$sites, rownames(scrs$sites), col = "black",
                 cex = 0.7, ...)
        }
        if(type == "points"){
            points(scrs$sites, col = "black", pch = 1, cex = 0.7, ...)
        }
    }
    if (axes) {
        axis(1)
        axis(2)
        box()
    }
    if(ann)
        title(main = main, sub = sub, xlab = xlab, ylab = ylab, ...)
    class(scrs) <- "ordiplot"
    invisible(scrs)
}
