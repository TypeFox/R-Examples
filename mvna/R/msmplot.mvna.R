msmplot <- function(object, ...) {
    UseMethod("msmplot")
}


###############################################
# Function to plot the multistate model
# along with the estimated cumulative hazards.
###

msmplot.mvna <- function(object, tr.choice, graph, layout = "dot",
                         recip.edges = "combined", order = "LR", unit = "npc",
                         width, height, just = "center",  
                         mgp = c(1.3, 0.4, 0), cex = 0.6, mtext = TRUE,
                         label.plots, side = 3, lcex = 0.9, conf.int = TRUE, 
                         level = 0.05, xlab = "Time", ylab = "", ...) {
    if (!inherits(object, "mvna")) {
        stop("Argument 'object' must be of class 'mvna'")
    }
    if (sum(c("Rgraphviz", "gridBase") %in% loadedNamespaces()) != 2) {
        stop("Packages 'Rgraphviz' and 'gridBase' must be installed and loaded")
    }
    if (!(layout %in% c("dot", "twopi", "neato", "circo", "fdp"))) {
        stop("Non existing layout method. See the Rgraphviz package help pages")
    }
    if (missing(tr.choice)) {
        tr.choice <- names(object)[1:(length(object) - 7)]
    }
    else {
        if (!all(tr.choice %in% names(object))) {
            stop("Argument 'tr.choice' and the possible transitions must match")
        }
    }
    x <- object[tr.choice]
    if (missing(width)) {
        width <- rep(0.25, length(x))
    }
    else if (length(width) < length(x)) {
        width <- width * rep(1, length(x))
    }
    if (missing(height)) {
        height <- rep(0.2, length(x))
    }
    else if (length(height) < length(x)) {
        height <- height * rep(1, length(x))
    }
    lim <- double(length(x))
    if (conf.int) {
        for (i in 1:length(x)) {
            x[[i]]$ciplus <- x[[i]]$na*exp((qnorm(1-level/2)*sqrt(x[[i]]$var1))/x[[i]]$na)
            x[[i]]$cimoins <- x[[i]]$na*exp((-qnorm(1-level/2)*sqrt(x[[i]]$var1))/x[[i]]$na)
            lim[i] <- max(x[[i]]$ciplus, na.rm = TRUE)
        }
    }
    else {
        lim <- sapply(1:length(x), function(i) {
            max(x[[i]]$ciplus, na.rm=TRUE)
        })
    }
    ylim <- c(0, max(lim, na.rm=TRUE))
    xlim <- c(0, as.numeric(dimnames(object$nev)[[3]][length(dimnames(object$nev)[[3]])]))
    oldpar <- par(no.readonly = TRUE)
### If missing graph, creation 
    if (missing(graph)) {
        msm <- ftM2graphNEL(as.matrix(object$trans))
        graph.msm <- plot(msm, layout, recipEdges = recip.edges,
                          attrs=list(graph=list(rankdir = order)))
    }
    else {
        plot(graph)
    }
    if (missing(label.plots)) label.plots <- tr.choice
### cumulative hazard plots
    vps <- baseViewports()
    pushViewport(vps$inner, vps$figure, vps$plot)
    cat("Order of the transitions: ", paste(tr.choice, "|"), "\n")
    for (i in 1:length(x)) {
        loc <- grid.locator(unit = unit)
        pushViewport(viewport(x = loc$x, y = loc$y,
                              width = unit(width[i], unit),
                              height = unit(height[i], unit),
                              just = just))
        par(plt = gridPLT(), new = TRUE, cex = cex, mgp = mgp)
        plot(x[[i]]$time, x[[i]]$na, type="s", xlim = xlim,
             ylim = ylim, xlab = xlab, ylab = ylab, ...)
        if (conf.int) {
            lines(x[[i]]$time, x[[i]]$ciplus, lty = 2, type= "s", ...)
            lines(x[[i]]$time, x[[i]]$cimoins, lty = 2, type = "s", ...)
        }
        if (mtext) {
            mtext(label.plots[i], side = side, cex = lcex)
        }
        popViewport()
    }
    popViewport(3)
    par(oldpar)
    object
}
