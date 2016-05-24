
###_ + Internal Function
".night" <- function(time, sunrise.time, sunset.time)
{
    ## Value: A list with sunset and sunrise times for dates in 'time'
    ## --------------------------------------------------------------------
    ## Arguments: Passed from plotTDR
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    morn.uniq <- unique(format(time, format=paste("%Y-%m-%d", sunrise.time)))
    tz <- ifelse(is.null(attr(time, "tzone")), "", attr(time, "tzone"))
    morn <- as.POSIXct(morn.uniq, tz=tz) + 86400
    morn.before <- morn[1] - 86400
    morn.all <- rbind(data.frame(x=morn.before), data.frame(x=morn))[[1]]
    night.uniq <- unique(format(time, format=paste("%Y-%m-%d", sunset.time)))
    night <- as.POSIXct(night.uniq, tz=tz)
    night.before <- night[1] - 86400
    night.all <- rbind(data.frame(x=night.before), data.frame(x=night))[[1]]
    list(sunrises=morn.all, sunsets=night.all)
}

###_ + Main Function
".plotTDR" <- function(time, depth, concurVars=NULL, xlim=NULL, depth.lim=NULL,
                       xlab="time (dd-mmm hh:mm)", ylab.depth="depth (m)",
                       concurVarTitles=deparse(substitute(concurVars)),
                       xlab.format="%d-%b %H:%M", sunrise.time="06:00:00",
                       sunset.time="18:00:00", night.col="gray60",
                       dry.time=NULL, phase.factor=NULL, plot.points=FALSE,
                       interact=TRUE, key=TRUE, cex.pts=0.4, ...)
{
    ## Value: Returns (invisibly) a list with coordinates for each zoc'ed
    ## time window.  Also Plot time, depth, and other concurrent data.
    ## --------------------------------------------------------------------
    ## Arguments: time=POSIXct; depth=numeric vector with depth readings,
    ## concurVars=matrix of numeric data with concurrent data to plot,
    ## xlim=POSIXct vector with lower and upper time limits to plot,
    ## depth.lim=vector with lower and upper depth limits, dry.time=subset
    ## of time corresponding to observations considered to be dry;
    ## phase.factor=factor classifying each reading, xlab=title for the x
    ## axis, ylab.depth=title for the depth axis, concurVarTitles=string
    ## vector with titles for the additional variables, xlab.format=format
    ## string for formatting time in x axis, sunrise.time=string specifying
    ## the time of sunrise, sunset.time=string specifying sunset time,
    ## night.col=color for masking night times, key=logical whether to draw
    ## a legend; ...=parameters passed to par.
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    nights <- .night(time, sunrise.time, sunset.time)
    nconcurVars <- ifelse(is.null(concurVars), 0, ncol(concurVars))
    plotrows <- nconcurVars + 1
    ncheight <- 1.35 * (1/plotrows)
    lheights <- c(rep((1 - ncheight)/nconcurVars, nconcurVars), ncheight)
    mardepthonly <- c(4, 4, 1, 1) + 0.1    # for depth plot only
    mardepthmore <- c(4, 4, -0.1, 1) + 0.1 # for depth plot and more vars
    martop <- c(-0.1, 4, 1, 1) + 0.1       # for top plot
    marnontop <- c(-0.1, 4, -0.1, 1) + 0.1 # for plots between depth and top one
    orig <- structure(0, class=class(time), tzone=attr(time, "tzone"))
    "plot.fun" <- function(xlim, ylim) {
        xticks <- orig + seq(from=xlim[1], to=xlim[2], length=20)
        if(is.null(concurVars)) {
            par(las=1, bty="n", mar=mardepthonly, ...)
        } else {
            par(las=1, bty="n", mar=mardepthmore, ...)
            layout(matrix(seq(plotrows, 1), nrow=plotrows, ncol=1),
                   heights=lheights)
        }
        now <- (time >= xlim[1]) & (time <= xlim[2])
        depth.now <- depth[now]
        time.now <- time[now]
        plot(depth.now ~ time.now, type="n", xlim=xlim, ylim=ylim,
             xlab=xlab, ylab=ylab.depth, xaxt="n", yaxt="n")
        usr <- par("usr")
        xleft <- pmax(unclass(nights$sunsets), usr[1])
        xright <- pmin(unclass(nights$sunrises), usr[2])
        rect(xleft, usr[3], xright, usr[4], col=night.col, border=NA)
        if (!is.null(dry.time)) segments(dry.time, usr[4], dry.time, usr[4],
                                         lwd=4, col="tan")
        axis.POSIXct(side=1, time.now, at=xticks, format=xlab.format)
        axis(side=2)
        lines(time.now, depth.now)
        if (!is.null(phase.factor)) {
            phase.factor <- phase.factor[now, drop=TRUE]
            nlevs <- nlevels(phase.factor)
            ncolors <- max(3, min(nlevs, 9))
            colors <- hsv(seq(0, 0.9, length=ncolors), 0.8, 0.95)
            points(time.now, depth.now, col=colors[phase.factor],
                   pch=19, cex=cex.pts)
            if (key && nlevs < 10 && nlevs > 0) {
                legend("bottomright", legend=levels(phase.factor), col=colors,
                       pch=19, cex=0.7, ncol=nlevs, bg="white")
            }
        } else if (plot.points) {
            points(time.now, depth.now, pch=19, cex=cex.pts)
        }
        if (!is.null(concurVars)) {
            if (length(concurVarTitles) != nconcurVars) {
                concurVarTitles <- rep(concurVarTitles, length.out=nconcurVars)
            }
            for (i in seq(nconcurVars)) {
                vari <- concurVars[now, i]
                if (i == nconcurVars) par(mar=martop) else par(mar=marnontop)
                ylim <- range(vari, na.rm=TRUE)
                plot(vari ~ time.now, type="n", xaxt="n", ylim=ylim,
                     xlab="", xlim=xlim, bty="n", ylab=concurVarTitles[i])
                usr <- par("usr")    # to watch out for change in y coords
                rect(xleft, usr[3], xright, usr[4], col=night.col, border=NA)
                lines(time.now, vari)
                if (!is.null(phase.factor)) { # we already have 'colors'
                    points(time.now, vari, col=colors[phase.factor], pch=19,
                           cex=cex.pts)
                } else if (plot.points) points(time.now, vari, pch=19, cex=cex.pts)
                axis(side=2)
            }
        }
    }
    if (!interact) {
        rx <- range(as.numeric(time))   # max and min of dates
        xlim <- if(is.null(xlim)) rx else as.numeric(xlim)
        ylim <- if (is.null(depth.lim)) {
            rev(range(depth, na.rm=TRUE))
        } else rev(depth.lim)
        plot.fun(xlim=xlim, ylim=ylim)
    } else {
        requireNamespace("tcltk", quietly=TRUE) ||
            stop("tcltk support is absent")
        rx <- range(as.numeric(time))   # max and min of dates
        diffrx <- diff(rx)
        xlim <- x10 <- if(is.null(xlim)) { # define xlim if not there already
            rx + (diffrx * 0.01)           # add 1% to each side
        } else as.numeric(xlim)
        xlmid <- xm0 <- mean(xlim)     # two vars with date range midpoint
        xr0 <- diff(xlim)                     # range of xlim
        xZoom <- tcltk::tclVar(100)           # initialize zoom factor
        xlmid <- tcltk::tclVar(xlmid)     # initialize date range midpoint
        xZ <- as.numeric(tcltk::tclvalue(xZoom)) # these 2 are to be dynamically changed
        xM <- as.numeric(tcltk::tclvalue(xlmid))
        ylim <- if (is.null(depth.lim)) {
            rev(range(depth, na.rm=TRUE)) * 1.1
        } else rev(depth.lim)
        yMax <- tcltk::tclVar(ylim[1])
        yTop <- as.numeric(tcltk::tclvalue(yMax))
        replot <- function(...) {
            xZ <<- as.numeric(tcltk::tclvalue(xZoom))
            xM <<- as.numeric(tcltk::tclvalue(xlmid))
            xr.half <- (xr0/2) * 100/xZ
            xlim <- xM + c(-xr.half, xr.half)
            yTop <<- as.numeric(tcltk::tclvalue(yMax))
            ylim <- c(yTop, ylim[2])
            plot.fun(xlim=xlim, ylim=ylim)
        }
        replot.maybe <- function(...) {
            if(as.numeric(tcltk::tclvalue(xZoom)) != xZ ||
               as.numeric(tcltk::tclvalue(xlmid)) != xM ||
               as.numeric(tcltk::tclvalue(yMax)) != yTop) replot()
        }
        coords <- list()
        zocrange <- function() {
            coords[[length(coords) + 1]] <<- locator(2)
            tcltk::tkgrab.release(base)
        }

        base <- tcltk::tktoplevel()
        tcltk::tkwm.title(base, "diveMove")
        tcltk::tkwm.deiconify(base)
        tcltk::tkgrab.set(base)
        tcltk::tkfocus(base)

        base.frame <- tcltk::tkframe(base, borderwidth=3)

        dep.frame <- tcltk::tkframe(base.frame, relief="groove", borderwidth=2)

        x.frame <- tcltk::tkframe(base.frame)
        xr.frame <- tcltk::tkframe(x.frame, relief="groove", borderwidth=2)
        xmid.frame <- tcltk::tkframe(x.frame, relief="groove", borderwidth=2)
        zoc.pts <- tcltk::tkbutton(base.frame, text="Zero-Offset\nCorrect a Range",
                                   command=zocrange)
        q.but <- tcltk::tkbutton(base.frame, text="Quit",
                                 command=function() tcltk::tkdestroy(base))

        ## Zoom
        diffx <- diff(as.numeric(time))
        diffxOK <- min(diffx[diffx > 0]) * 40 # zoom up to 40 observations
        maxZoom <- (diffrx / diffxOK) * 100 # maximum zoom depends on time range
        tzoom.l <- tcltk::tklabel(xr.frame, text="Date Zoom (%)")
        tzoom.s <- tcltk::tkscale(xr.frame, command=replot.maybe, from=100,
                                  to=maxZoom, showvalue=TRUE, variable=xZoom,
                                  resolution=100, length=200, orient="horiz")
        ## Pan
        tpan.l <- tcltk::tklabel(xmid.frame, text="Pan through Date")
        tpan.s <- tcltk::tkscale(xmid.frame, command=replot.maybe,
                                 from=xm0 - xr0, to=xm0 + xr0,
                                 showvalue=FALSE, variable=xlmid,
                                 resolution=xr0/2000, length=200,
                                 orient="horiz")
        ## Maximum depth selection
        maxdep.l <- tcltk::tklabel(dep.frame, text="Max. Depth (m)")
        maxdep.s <- tcltk::tkscale(dep.frame, command=replot.maybe,
                                   from=0, to=ylim[1], length=150,
                                   showvalue=TRUE, variable=yMax,
                                   orient="vertical")

        ## Grid all the widgets together
        tcltk::tkgrid(base.frame)
        tcltk::tkgrid(dep.frame, rowspan=3, column=0)
        tcltk::tkgrid(maxdep.l); tcltk::tkgrid(maxdep.s, sticky="ns")
        tcltk::tkgrid(x.frame, row=0, column=1, columnspan=2, sticky="n")
        tcltk::tkgrid(xr.frame)
        tcltk::tkgrid(tzoom.l, sticky="ew"); tcltk::tkgrid(tzoom.s, sticky="ew")
        tcltk::tkgrid(xmid.frame)
        tcltk::tkgrid(tpan.l, sticky="ew"); tcltk::tkgrid(tpan.s, sticky="ew")
        tcltk::tkgrid(zoc.pts, row=2, column=1, sticky="ns")
        tcltk::tkgrid(q.but, row=2, column=2)

        if (getRversion() >= "2.14.2") replot()
        tcltk::tkwait.window(base)
        invisible(coords)
    }
}



###_ + Emacs local variables
## Local variables:
## allout-layout: (+ : 0)
## End:
