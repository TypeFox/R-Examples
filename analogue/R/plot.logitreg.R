`plot.logitreg` <- function(x, group = "all",
                            npred = 100,
                            conf.int = 0.9,
                            conf.type = c("none","polygon","lines"),
                            xlab = expression(D[ij]),
                            ylab = "Pr (A+ | d)",
                            rug = TRUE,
                            ticksize = 0.02,
                            col = "red",
                            ref.col = "lightgrey",
                            lwd = 2,
                            conf.lwd = 1,
                            conf.lty = "dashed",
                            shade = "lightgrey",
                            ...) {
    pfun <- function(x, min, max, npred, conf.int, conf.type, rug,
                     col, lwd, shade, main, ticksize, ref.col,
                     conf.lwd, conf.lty, ...) {
        dat <- data.frame(Dij = seq(min, max, length = npred))
        ilogit <- family(x)$linkinv
        pred <- predict(x, newdata = dat, type = "link", se.fit = TRUE)
        dat$pred <- ilogit(pred$fit)
        crit.t <- qt(conf.int - ((1- conf.int) / 2), df = npred - 2)
        dat$upper <- ilogit(pred$fit + (crit.t * pred$se.fit))
        dat$lower <- ilogit(pred$fit - (crit.t * pred$se.fit))
        ## add a bit of extra space
        y.lim <- c(-0.02, 1.02)
        plot(x = dat$Dij, y = dat$pred, type = "n",
             ylim = y.lim,ylab = "", xlab = "",
             axes = FALSE, main = main, ...)
        ## draw reference line
        abline(h = c(0,1), col = ref.col)
        ## draw rug plots?
        if(rug) {
            rug(x$data$Dij[x$data$analogs], side = 3, ticksize = ticksize)
            rug(x$data$Dij[!x$data$analogs], side = 1, ticksize = ticksize)
        }
        ## draw confidence interval?
        if(conf.type == "polygon") {
            polygon(x = c(dat$Dij, rev(dat$Dij)),
                    y = c(dat$upper, rev(dat$lower)),
                    col = shade, border = shade)
        } else if(conf.type == "lines") {
            lines(x = dat$Dij, y = dat$upper, col = col, lwd = conf.lwd,
                  lty = conf.lty)
            lines(x = dat$Dij, y = dat$lower, col = col, lwd = conf.lwd,
                  lty = conf.lty)
        }
        ## draw fitted function
        lines(x = dat$Dij, y = dat$pred, col = col, lwd = lwd)
        axis(1)
        axis(2)
        box()
    }
    ## process group; allow 'all' or one of the actual groups
    if (length(group) > 1) {
        group <- group[1]
        warning("More than 1 'group' specified. Using only the first\nDid you mean to use '\"all'\"")
    }
    n.groups <- length(x$models)
    g.names <- names(x$models)
    GROUP <- c("all", g.names)
    group <- match.arg(group, GROUP)
    conf.type <- match.arg(conf.type)
    if (group == "all") {
        n.plot <- n.groups
        xy.nums <- n2mfrow(n.plot)
        layout(matrix(seq_len(prod(xy.nums)), nrow = xy.nums[1],
            ncol = xy.nums[2]))
        op <- par(mar = c(2, 2, 3, 1) + 0.1,
                  oma = c(2, 2, 2, 0), no.readonly = TRUE)
        on.exit({
            par(op)
            layout(1)
        })
        min <- 0
        max <- max(x$models[["Combined"]]$data$Dij)
    }
    else {
        n.plot <- 1
        min <- 0
        max <- max(x$models[[group]]$data$Dij)
        ##xy.nums <- rep(1, 2)
    }
    for (i in seq_len(n.groups)) {
        if (group != "all" && group != g.names[i])
            next
        pfun(x$models[[i]], min = min, max = max, npred = npred,
             col = col, lwd = lwd, shade = shade, main = g.names[i],
             conf.type = conf.type, conf.int = conf.int, rug = rug,
             ticksize = ticksize, ref.col = ref.col,
             conf.lty = conf.lty, conf.lwd = conf.lwd, ...)
    }
    if (n.plot > 1) {
        title(xlab = xlab, outer = TRUE, line = 0.5, cex.lab = 1.3)
        title(ylab = ylab, outer = TRUE, line = 0.5, cex.lab = 1.3)
        title(main = "Posterior probability of analogue", outer = TRUE,
            line = 0.5, cex.main = 1.3)
    }
    else {
        title(xlab = xlab, ylab = ylab, sub = "")
    }
    invisible(x)
}
