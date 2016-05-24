plot.etmCIF <- function(x, which.cif, xlim, ylim,
                        ylab = "Cumulative Incidence", xlab = "Time",
                        col = 1, lty, lwd = 1, ci.type = c("none", "bars", "pointwise"),
                        ci.fun = "cloglog", ci.col = col, ci.lty = 3,
                        legend = TRUE, legend.pos, curvlab, legend.bty = "n",
                        pos.ci = 27, ci.lwd = 3, 
                        ...) {

    if (!inherits(x, "etmCIF")) {
        stop("'x' must be of class 'etmCIF'")
    }
    ci.type <- match.arg(ci.type)
    
    tr.choice <- paste(x[[1]]$trans[, 1], x[[1]]$trans[, 2])
    l.x <- NCOL(x$X)
    n.trans <- length(tr.choice)

    if (missing(which.cif)) {
        tr.choice <- paste(0, x$failcode, sep = " ")
    } else {
        tr.choice <- paste(0, which.cif, sep = " ")
        ## A small test on tr.choice
        ref <- sapply(1:length(x[[1]]$state.names), function(i) {
            paste(x[[1]]$state.names, x[[1]]$state.names[i])
        })
        ref <- matrix(ref)
        if (sum(tr.choice %in% ref == FALSE) > 0)
            stop("Argument 'which.cif' and causes of failure must match")
    }

    n.what <- length(tr.choice)
    
    max.time <- max(sapply(x[1:l.x], function(ll) {
        max(ll$time)
    }))

    if (missing(ylim)) ylim <- c(0, 1)
    if (missing(xlim)) xlim <- c(0, max.time)
    if (missing(lty)) {
        lty <- seq_len(n.what * l.x)
    }
    else if (length(lty) < (l.x * n.what)) {
        lty <- lty * rep(1, l.x * n.what)
    }
    if (length(col) < l.x * n.what)
        col <- col * rep(1, l.x * n.what)

    conf.int <- if (ci.type == "pointwise") TRUE else FALSE
    if (ci.type != "none") {
        if (missing(ci.col)) {
            ci.col <- col
        } else {
            if (length(ci.col) < (l.x * n.what)) {
                ci.col <- ci.col * rep(1, l.x * n.what)
            }
        }
        if (missing(ci.lty)) {
            ci.lty <- lty
        } else {
            if (length(ci.lty) < (l.x * n.what)) {
                ci.lty <- ci.lty * rep(1, l.x * n.what)
            }
        }
    }

    plot(xlim, ylim, xlab = xlab, ylab = ylab,
         xlim = xlim, ylim = ylim, type = "n", ...)
    
    summx <- lapply(x[1:l.x], summary, ci.fun = ci.fun)

    if (length(pos.ci) < l.x)  pos.ci <- rep(pos.ci, l.x)

    for (i in seq_len(l.x)) {
        for (j in seq_along(tr.choice)) {
            lines(x[[i]], tr.choice = tr.choice[j],
                  col = col[j + (i - 1) * n.what], lty = lty[j + (i - 1) * n.what],
                  lwd = lwd, conf.int = conf.int,...)

            if (ci.type == "bars") {
                ind <- findInterval(pos.ci[i], summx[[i]][[tr.choice[j]]]$time)
                segments(pos.ci[i], summx[[i]][[tr.choice[j]]]$lower[ind],
                         pos.ci[i], summx[[i]][[tr.choice[j]]]$upper[ind], 
                         lwd = ci.lwd, col = ci.col[j + (i - 1) * n.what],
                         lty = ci.lty[j + (i - 1) * n.what],...)
            }
        }
    }
        
    if (legend) {
        if (missing(legend.pos)) {
            legend.pos <- "topleft"
        }
        if (missing(curvlab)) {
            cdc <- sapply(strsplit(sub("\\s", "|", tr.choice), "\\|"),
                                              "[", 2)
##            cdc <- sapply(strsplit(tr.choice, " "), "[", 2)
            if (l.x == 1) {
                curvlab <- paste("CIF ", cdc, sep = "")
            } else {
                if (length(cdc) == 1) {
                    curvlab <- paste("CIF ", cdc, "; ", rownames(x$X), "=", x$X, sep = "")
                } else {
                    curvlab <- as.vector(sapply(seq_along(x$X), function(j){
                        paste("CIF ", cdc, "; ", rownames(x$X), "=", x$X[j], sep = "")
                    }))
                }
            }
        }
        if (is.list(legend.pos)) legend.pos <- unlist(legend.pos)
        if (length(legend.pos) == 1) {
            xx <- legend.pos
            yy <- NULL
        }
        if (length(legend.pos) == 2) {
            xx <- legend.pos[1]
            yy <- legend.pos[2]
        }
        args <- list(...)
        ii <- pmatch(names(args),
                     names(formals("legend")[-charmatch("bty",names(formals("legend")))]))
        do.call("legend", c(list(xx, yy, curvlab, col=col, lty=lty, lwd = lwd, bty = legend.bty),
                            args[!is.na(ii)]))
    }
    
    invisible()
    
}
