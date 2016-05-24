plot.roc <- function(x,
                     which = c(1:3,5),
                     group = "Combined",
                     prior = NULL,
                     show.stats = TRUE,
                     abline.col = "grey",
                     abline.lty = "dashed",
                     inGroup.col = "red",
                     outGroup.col = "blue",
                     lty = "solid",
                     caption = c("ROC curve",
                       "Dissimilarity profiles",
                     "TPF - FPF vs Dissimilarity",
                     "Likelihood ratios"),
                     legend = "topright",
                     ask = prod(par("mfcol")) < length(which) &&
                     dev.interactive(),
                     ...) {
    if(!inherits(x, "roc"))
        stop("Plot method only for objects of class \"roc\".")
    if (!is.numeric(which) || any(which < 1) || any(which > 5))
        stop("'which' must be in 1:5")
    show <- rep(FALSE, 5)
    show[which] <- TRUE
    one.fig <- prod(par("mfcol")) == 1
    if(ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    if(any(show[4:5])){
        l.ratios <- bayesF(x, prior = prior)
        if(group != "Combined") {
            g.names <- names(l.ratios[seq_len(length(l.ratios) - 2)])
            group <- match.arg(group, g.names)
        }
    }
    method <- attr(x, "method")
    ## replicate lty if length 1
    if(length(lty) == 1L)
      lty <- rep(lty, 2)
    if(show[1]) {
        plot(x$roc[[group]]$FPE, x$roc[[group]]$TPF, type = "n",
             ylab = "TPF (sensitivity)",
             xlab = "1 - TNF (1 - specificity)")
        abline(0, 1, col = abline.col)
        lines(x$roc[[group]]$FPE, x$roc[[group]]$TPF, lty = lty[1L], ...)
        mtext(caption[1], side = 3, line = 1.7, font = 2)
        if(show.stats) {
            txt <- paste("AUC =", round(x$roc[[group]]$AUC, 3))
            legend("bottomright", legend = txt, pch = NA,
                   bty = "n", cex = 0.8)
        }
        box()
    }
    if(show[2]) {
        dens.in <- density(x$roc[[group]]$analogue$yes, from = 0)
        dens.out <- density(x$roc[[group]]$analogue$no, from = 0)
        xlims <- switch(method,
                        SQchord = c(0,2),
                        chord = c(0, sqrt(2)),
                        c(0, max(x$dissims)))
        ylims <- range(0, dens.in$y, dens.out$y)
        plot(dens.in$x, dens.in$y, type = "n", axes = FALSE,
             xlim = xlims, ylim = ylims,
             ylab = "Density",
             xlab = paste("Dissimilarity (", attr(x, "method"), ")",
             sep = ""))
        abline(h = 0, col = abline.col)
        lines(dens.in, col = inGroup.col, lty = lty[1L], ...)
        lines(dens.out, col = outGroup.col, lty = lty[2L], ...)
        abline(v = x$roc[[group]]$optimal, lty = abline.lty,
               col = abline.col)
        axis(side = 2)
        axis(side = 1)
        box()
        mtext(caption[2], side = 3, line = 1.7, font = 2)
        legend("topright", legend = c("Analogue", "Not Analogue"),
               col = c(inGroup.col, outGroup.col),
               lty = lty[c(1L,2L)],
               bty = "n", cex = 0.8,
               inset = 0.01, ...)
    }
    if(show[3]) {
        cutpoints <- unique(rev(x$roc[[group]]$roc.points))
        roc.values <- x$roc[[group]]$TPF - x$roc[[group]]$FPE
        plot(cutpoints, roc.values, type = "n",
             ylab = "TPF - (1 - TNF)",
             xlab = paste("Dissimilarity (", attr(x, "method"), ")",
             sep = ""))
        abline(h = 0, col = abline.col)
        lines(cutpoints, roc.values, lty = lty[1L], ...)
        abline(v = x$roc[[group]]$optimal, lty = abline.lty,
               col = abline.col)
        mtext(caption[3], side = 3, line = 1.7, font = 2)
        box()
    }
    if(show[4]) {
        pos <- l.ratios[[group]]$bayesF$pos
        neg <- l.ratios[[group]]$bayesF$neg
        dissims <- rev(x$roc[[group]]$roc.points)
        plot(dissims, pos, type = "n", axes = FALSE, ylab = "LR (+)",
             xlab = paste("Dissimilarity (", attr(x, "method"), ")",
             sep = ""))
        abline(v = x$roc[[group]]$optimal, lty = abline.lty,
               col = abline.col)
        lines(dissims, pos, col = inGroup.col, lty = lty[1L], ...)
        axis(side = 1)
        axis(side = 2)
        usr <- par("usr")
        finite.vals <- is.finite(neg)
        rany <- (max(neg[finite.vals]) -
                 min(neg[finite.vals])) * 0.04
        par(usr = c(usr[1:2], min(neg[finite.vals]) - rany,
            max(neg[finite.vals]) + rany))
        lines(dissims, neg, col = outGroup.col, lty = lty[2L], ...)
        axis(side = 4)
        box()
        mtext(caption[4], side = 3, line = 1.7, font = 2)
        legend("topright", legend = c("LR (+)", "LR (-)"),
               col = c(inGroup.col, outGroup.col),
               lty = lty,
               bty = "n", cex = 0.8,
               inset = 0.01, ...)
    }
    if(show[5]) {
        plot(l.ratios, group = group,
             abline.col = abline.col, col = inGroup.col,
             ylab = "Pr (A+ | d)",
             xlab = paste("Dissimilarity (", attr(x, "method"), ")",
             sep = ""), lty = lty[1L], ...)
        ##mtext(caption[5], side = 3, line = 1.7, font = 2)
    }
    invisible()
}
