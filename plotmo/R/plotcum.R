# plotcum.R

plotmo_cum <- function(rinfo, info, nfigs=1, add=FALSE,
                       cum.col1, grid.col, jitter=0, cum.grid="percentages", ...)
{
    trans.resids <- abs(rinfo$scale * rinfo$resids)
    # TODO what happens here if NA in trans.resids (leverage==1)
    ecdf <- ecdf(trans.resids[,1])
    xlab <- rinfo$name
    xlab <- sprintf("abs(%ss)", xlab)
    cum.grid <- match.choices(cum.grid, c("none", "grid", "percentages"))
    annotation.cex <- .7 * dot("cum.cex", DEF=1, ...)
    if(!add && info && cum.grid == "percentages") {
        # ensure right margin big enough for right hand labels
        old.mar <- par("mar")
        if(old.mar[4] < 3.5) {
            on.exit(par(mar=old.mar))
            par(mar=c(old.mar[1:3], annotation.cex * 5))
        }
    }
    if(is.na(cum.col1))
        cum.col1 <- dot("cum.col", DEF=1, ...)
    cum.col1 <- cum.col1[1] # no recycling

    # user can set xlim and ylim if this is the only figure
    xlim <- dot("xlim", DEF=NULL, ...)
    if(nfigs > 1 || !is.specified(xlim))
        xlim <- range(abs(rinfo$scale * rinfo$resids), na.rm=TRUE)
    xlim <- fix.lim(xlim)
    ylim <- dot("ylim", DEF=NULL, ...)
    if(nfigs > 1 || !is.specified(ylim))
        ylim <- c(ylim=if(info) -.1 else 0,
                  ymax=if(cum.grid == "percentages") 1 + annotation.cex * .06
                       else                          1)
    ylim <- fix.lim(ylim)

    call.plot(stats::plot.stepfun, PREFIX="cum.", drop.cum.grid=1,
        force.x          = ecdf,
        force.add        = add,
        force.main       = dot("main", DEF="Cumulative Distribution", ...),
        force.xlim       = xlim,
        force.ylim       = ylim,
        force.xlab       = xlab,
        force.ylab       = "Proportion",
        force.col.points = NA, # finer resol graph (points are big regardless of pch)
        force.col        = cum.col1,
        force.col.hor    = cum.col1,
        force.col.vert   = cum.col1,
        ...)
    if(!add) {
        if(info)
            draw.density.along.the.bottom(abs(trans.resids), ...)
        if(cum.grid %in% c("grid", "percentages")) {
            linecol <- if(is.specified(grid.col)) grid.col else "lightgray"
            # add annotated grid lines, unattractive but useful
            for(h in c(0,.25,.5,.75,.90,.95,1)) # horizontal lines
                abline(h=h, lty=1, col=linecol)
            probs <- c(0, .25, .50, .75, .9, .95, 1)
            q <- quantile(trans.resids,
                          probs=probs, na.rm=TRUE)
            for(v in q)    # vertical lines at 0,25,50,75,90,95,100% quantiles
                abline(v=v, lty=1, col=linecol)
            box() # abline overwrite the box, so restore it
            if(cum.grid == "percentages") {
                draw.percents.on.top(probs, q, annotation.cex)
                if(info)
                    draw.quantiles.on.right.side(probs, q, annotation.cex)
            }
            # replot data over grid
            call.plot(stats::plot.stepfun, PREFIX="cum.", drop.cum.grid=1,
                force.x          = ecdf,
                force.add        = TRUE,
                force.xlim       = xlim,
                force.col.points = NA,
                force.col        = cum.col1,
                force.col.hor    = cum.col1,
                force.col.vert   = cum.col1,
                ...)
        }
    }
}
# Adding percents and quantiles on the wrong axes is considered a no no,
# but here we are more-or-less forced into it because the percentile text
# can be too long to display on the "correct" axis.
draw.percents.on.top <- function(probs, q, annotation.cex)
{
    is.space.available <- function(i) # is horizontal space available
    {
        q[i] - q[i-1] > 1.2 * strwidth && q[i+1] - q[i] > 1.2 * strwidth
    }
    draw.percent <- function(i, label)
    {
        # xpd=NA to allow text out of plot region (usually not needed)
        x <- q[i]
        if(i == 1)
            x <- x + .05 * strwidth # so 0% doesn't overwrite box
        else if(i == 7)
            x <- x - .3 * strwidth  # so 100% doesn't overwrite box
        text.on.white(x=x, y=1.05, label, annotation.cex, xmar=0, xpd=NA)
    }
    #--- draw.percents starts here ---
    strwidth <- strwidth("25%", cex=annotation.cex)
    draw.percent(1, "0%")
    if(is.space.available(2)) draw.percent(2, "25%")
    draw.percent(3, "50%")
    if(is.space.available(4)) draw.percent(4, "75%")
    draw.percent(5, "90%")
    if(is.space.available(6)) draw.percent(6, "95%")
    draw.percent(7, "100%")
}
draw.quantiles.on.right.side <- function(probs, q, annotation.cex)
{
    y <- TeachingDemos::spread.labs(x=probs,
            mindiff=1.2 * annotation.cex * strheight("A"), min=-.1)
    q[q < max(q) / 1e4] <- 0 # prevent labels like 2.22e-16
    text(1.01 * par("usr")[2], y, sprintf("%.3g", q),
         xpd=TRUE, cex=annotation.cex, adj=0)
}
