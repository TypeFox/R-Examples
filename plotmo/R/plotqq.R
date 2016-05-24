# plotqq.R

plotmo_qq <- function(rinfo, info, nfigs,
                      grid.col, smooth.col,
                      id.n, iresids, npoints, ...)
{
    old.pty <- par("pty")
    par(pty="s") # square
    on.exit(par(pty=old.pty))
    # we figure out the shape of the qq line with all resids but
    # plot only npoints points (selecting them with iresids)
    resids <- rinfo$scale * rinfo$resids
    # qqnorm sets NAs in trans.resids (leverage==1) to NA in
    # qq$x and qq$y, and thus NAs don't get plotted (R PR#3750)
    main <- dot("main", DEF=sprintf("%s QQ", rinfo$name), ...)
    qq <- qqnorm(resids, main=main, plot.it=FALSE)
    id.indices <- get.id.indices(resids, id.n)
    xlim <- NULL
    ylim <- NULL
    if(nfigs == 1) { # user can set xlim and ylim if this is the only figure
        xlim <- dot("xlim", DEF=NULL, ...)
        ylim <- dot("ylim", DEF=NULL, ...)
    }
    if(!is.specified(xlim) && !is.null(id.indices)) { # extra space for point labs?
        min <- min(qq$x, na.rm=TRUE)
        max <- max(qq$x, na.rm=TRUE)
        xlim <- c(min - .1 * (max - min), max + .1 * (max - min))
    }
    if(!is.specified(ylim)) {
        min <- min(qq$y, na.rm=TRUE)
        max <- max(qq$y, na.rm=TRUE)
        ylim <- c(min, max)
        if(!is.null(id.indices))    # extra space for point labs?
            ylim <- c(min - .05 * (max - min), max + .05 * (max - min))
        if(info)                    # extra space for density plot?
            ylim[1] <- ylim[1] - .1 * (max - min)
    }
    xlim <- fix.lim(xlim)
    ylim <- fix.lim(ylim)

    # allow col.response as an argname for compat with old plotmo
    pt.col <- dot("col.response col.resp", DEF=1, ...)
    pt.col <- dot("pt.col col.points col.point col.residuals col.resid col",
                     EX=c(0,1,1,1,1,1), DEF=pt.col, NEW=1, ...)
    pt.col <- dot("qq.col col.residuals col.resid col", EX=c(0,1,1,1), DEF=pt.col, NEW=1, ...)
    # recycle
    pt.col <- repl(pt.col, length(resids))

    pt.cex <- dot("response.cex cex.response", DEF=1, ...)
    pt.cex <- dot("pt.cex cex.points cex.point", EX=c(0,1,1), DEF=pt.cex, NEW=1, ...)
    pt.cex <- dot("qq.cex cex.qq cex.residuals cex", EX=c(0,1,1,1), DEF=pt.cex, NEW=1, ...)
    pt.cex <- pt.cex * pt.cex(length(resids), npoints)
    pt.cex <- repl(pt.cex, length(resids))

    pt.pch <- dot("response.pch pch.response", DEF=20, ...)
    pt.pch <- dot(
        "qq.pch pt.pch pch.points pch.point pch.residuals pch",
        EX=c(1,0,0,1,1,1), DEF=pt.pch, NEW=1, ...)
    pt.pch <- repl(pt.pch, length(resids))

    ylab <- rinfo$name
    ylab <- sprintf("%s Quantiles", ylab)
    drop.line.col <- function(..., qqline.col=NA, qqline.lwd=NA, qqline.lty=NA)
    {
        call.plot(graphics::plot, PREFIX="qq.",
                  force.x    = qq$x[iresids],
                  force.y    = qq$y[iresids],
                  force.col  = pt.col[iresids],
                  force.cex  = pt.cex[iresids],
                  force.pch  = pt.pch[iresids],
                  force.main = main,
                  force.xlab = "Normal Quantiles",
                  force.ylab = ylab,
                  force.xlim = xlim,
                  force.ylim = ylim,
                  ...)
    }
    drop.line.col(...)
    if(is.specified(grid.col))
        grid(col=grid.col, lty=1)
    qqline.col <- dot("qqline.col", DEF="gray", ...)
    qqline.lwd <- dot("qqline.lwd", DEF=1,      ...)
    qqline.lty <- dot("qqline.lty", DEF=1,      ...)
    if(is.specified(qqline.col) &&
       is.specified(qqline.lwd) &&
       is.specified(qqline.lty))
        call.plot(qqline, force.y=resids,
                  force.col=qqline.col, force.lwd=qqline.lwd, force.lty=qqline.lty, ...)
    if(info) {
        # draw actual and theoretical density along the bottom
        usr <- par("usr") # xmin, xmax, ymin, ymax
        scale <- .1 * (usr[4] - usr[3]) / (max(qq$y) - min(qq$y))
        draw.density.along.the.bottom(qq$x, den.col=smooth.col, scale=scale, ...)
        draw.density.along.the.bottom(
            resids / sd(resids, na.rm=TRUE), # TODO correct?
            scale=scale, ...)
        legend("bottomright", inset=c(0,.06),
               legend=c("actual", "normal"),
               cex=.8, lty=1, col=c("gray57", smooth.col),
               box.col="white", bg="white", x.intersp=.2, seg.len=1.5)
    }
    if(is.specified(grid.col) || is.specified(qqline.col) || info) {
        # replot box and points because they may have been obscured
        box()
        drop.line.col <- function(..., qqline.col=NA, qqline.lwd=NA, qqline.lty=NA)
        {
            call.plot(graphics::points, PREFIX="qq.",
                      force.x    = qq$x[iresids],
                      force.y    = qq$y[iresids],
                      force.col  = pt.col[iresids],
                      force.cex  = pt.cex[iresids],
                      force.pch  = pt.pch[iresids],
                      ...)
        }
        drop.line.col()
    }
    if(!is.null(id.indices))
        plotrix::thigmophobe.labels(
            x      = qq$x[id.indices], y=qq$y[id.indices],
            labels = rinfo$labs[id.indices],
            offset = .33, xpd=NA,
            font   = dot("label.font", DEF=1, ...)[1],
            cex    = .8 * dot("label.cex", DEF=1, ...)[1],
            col    = dot("label.col",
                        DEF=if(is.specified(smooth.col)) smooth.col else 2, ...)[1])
}
