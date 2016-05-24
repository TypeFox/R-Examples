"panelplot" <-
function (data, panel = points, totrows = 3, totcols = 2, oma = rep(2.5,
                                                              4),
par.strip.text = NULL)
{
    opar <- par(mfrow = c(totrows, totcols), mar = rep(0, 4),
                oma = oma)
    on.exit(par(opar))
    if (!is.null(par.strip.text)) {
        cex.strip <- par.strip.text$cex
        stripfac <- (par()$cin[2]/par()$pin[2]) * cex.strip *
            1
    }
    else stripfac <- 0
    fac <- names(data)
    if (is.null(fac))
        fac <- 1:length(data)
    nseq <- 1:length(fac)
    plot.new()
    for (index in nseq) {
        ilev <- fac[index]
        lis <- data[[ilev]]
        i <- totrows - ((index - 1)%/%totcols)
        j <- (index - 1)%%totcols + 1
        par(mfg = c(i, j, totrows, totcols))
        xlim <- lis$xlim
        ylim <- lis$ylim
        if (stripfac > 0) {
            strip.text <- fac[index]
            ylim[2] <- ylim[2] + diff(ylim) * stripfac
        }
        else strip.text <- NULL
#        plot.new()
        plot.window(unique(xlim), unique(ylim))
        if (!is.null(strip.text)) {
            chh <- par()$cxy[2]
            ht <- par()$usr[4] - 0.725 * chh
            abline(h = ht)
            xmid <- mean(par()$usr[1:2])
            text(xmid, ht + chh * 0.35, strip.text, cex = cex.strip)
        }
        panel(lis, nrows = i, ncols = j)
        box()
    }
}
