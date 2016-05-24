
RCplot <-
    function(x, r2 = residuals(x, "deviance")^2,
             alphabet = x$alpha, lab.horiz = k <= 20, do.call = TRUE,
             cex.axis = if(k <= 20) 1 else if(k <= 40) 0.8 else 0.6,
             y.fact = if(.Device == "postscript") 1.2 else 0.75,
             col = "gray70", xlab = "Context", main = NULL,
             med.pars = list(col = "red", pch = 12, cex= 1.25 * cex.axis),
             ylim = range(0, r2, finite=TRUE), ...)
{
    ## Author: Martin Maechler, Date:  1 Mar 2002, 17:39
    namx <- deparse(substitute(x))
    if(!is.vlmc(x)) stop("'x' must be a fitted VLMC object")
    fID <- id2ctxt(predict(x, type="id"), alpha = alphabet)
    ok <- fID != "NA"
    ## drop those with "NA" context (at least the first one!)
    ## FIXME: should we tell about this ?
    fID <- as.factor(fID[ok])
    r2 <- r2[ok]
    tfID <- table(fID)
    k <- length(tfID)
    if(is.null(main))
        main <- paste(if(!do.call)"VLMC", "Residuals vs. Context")
    labs <- c("#{obs}:", tfID)
    if(!lab.horiz && missing(ylim)) { ## use space *below* 0 line
        ## find out about string width: -> use x-direction
        op <- par(xaxs = "i"); plot.new(); plot.window(xlim=ylim, ylim=0:1)
        y0 <- y.fact * max(strwidth(labs, cex=cex.axis))
        ylim[1] <- min(ylim[1], - 1.4 * y0)
        par(op)
    }
    op <- par(cex.axis = cex.axis)
    on.exit(par(op))
    ## plot.factor calling (and returning) boxplot():
    rp <- plot(fID, r2, varwidth = TRUE, xlab = xlab, main = main,
               ylab = paste("residuals(",namx,", \"deviance\") ^ 2", sep=""),
               ylim = ylim, col = col, las = if(lab.horiz) 0 else 2, ...)
    abline(h = 0, lty = 3)
    if(any(i0 <- abs(meds <- rp$stats[3,]) < 1e-3))
        do.call("points", c(list(x=which(i0), y= meds[i0]), med.pars))
    if(lab.horiz)
        text(c(.2, 1:k), -.1, labs, xpd=FALSE, cex=cex.axis)
    else {
        text(c((par("usr")[1]+1)/2, 1:k), -0.3*y0, labs,
             xpd=NA, srt = 90, adj = c(1, 0.5), cex=cex.axis)
    }
    if(do.call) mtext(deparse(x$call))
    invisible(list(k = k, fID = fID, rp = rp))
}
