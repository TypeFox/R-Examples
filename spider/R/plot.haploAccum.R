plot.haploAccum <-
    function(x, add = FALSE, ci = 2, ci.type = c("bar","line","polygon"), 
             col = par("fg"), ci.col = col, ci.lty = 1, xlab,
             ylab = "Haplotypes", ylim, main = paste(x$method, "method of haplotype accumulation", sep=" "), ...)
{
    xaxvar <- x[["sequences"]]
    if (missing(xlab)) xlab <- "Sequences"
    ci.type <- match.arg(ci.type)
    if (!add) {
        if (missing(ylim))
            ylim <- c(1, max(x$n.haplotypes, x$n.haplotypes + ci*x$sd))
        plot(xaxvar, x$n.haplotypes, xlab=xlab, ylab=ylab, ylim=ylim,
             type="n", main = main, ...)
    }
    if (!is.null(x$sd) && ci)
        switch(ci.type,
               bar = segments(xaxvar, x$n.haplotypes - ci*x$sd, xaxvar,
                  x$n.haplotypes + ci*x$sd, col=ci.col, lty=ci.lty, ...),
               line = matlines(xaxvar, x$n.haplotypes + t(rbind(-ci,ci) %*% x$sd),
                 col=ci.col, lty=ci.lty, ...),
               polygon = polygon(c(xaxvar, rev(xaxvar)),
                 c(x$n.haplotypes - ci*x$sd, rev(x$n.haplotypes + ci*x$sd)), col=ci.col,
                 lty=ci.lty, main = main,  ...)
               )
    lines(xaxvar, x$n.haplotypes,col=col, ...)
    invisible()
}

