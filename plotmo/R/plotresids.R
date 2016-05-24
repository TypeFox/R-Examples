# plotresids.R

plotresids <- function(
    object,
    which,
    info,
    standardize,
    level,
    versus1,

    id.n,
    smooth.col,
    grid.col,
    jitter,

    npoints,
    center,

    type,

    fitted,
    rinfo,
    rsq,
    iresids,
    nversus,
    colname.versus1,
    force.auto.resids.xylim,

    SHOWCALL=NA, # this is here to absorb SHOWCALL from dots

    ...)
{
    stopifnot(length(which) == 1)
    info <- check.boolean(info)

    ok <- which %in% c(W3RESID,W5ABS:W9LOGLOG)
    if(!all(ok))
        stop0("which=", which[!ok][1], " is not allowed")

    id.indices <- NULL
    if(which %in% c(W3RESID, W4QQ:W8CUBE))
        id.indices <-
            get.id.indices(rinfo$scale * rinfo$resids, id.n,
                           if(nversus == V4LEVER)
                               hatvalues1(object, sprintf("versus=%g", V4LEVER))
                           else
                               NULL)

    level <- check.level.arg(level, zero.ok=TRUE)
    if(which %in% (W5ABS:W9LOGLOG))
        level <- 0 # no pints

    pints <- NULL
    cints <- NULL
    level.shade  <- dot("level.shade  shade.pints", DEF="mistyrose2", ...)
    level.shade2 <- dot("level.shade2 shade.cints", DEF="mistyrose4", ...)
    if(which == W3RESID && is.specified(level)) {
        p <- plotmo.pint(object, newdata=NULL, type, level, trace=0)
        stopifnot(is.null(p$fit) || (p$fit - fitted == 0))
        if(is.specified(level.shade) && !is.null(p$upr)) {
            pints <- data.frame(upr=rinfo$scale * (p$upr - fitted),
                                lwr=rinfo$scale * (p$lwr - fitted))
            colnames(pints) <- c("upr", "lwr")
        }
        if(is.specified(level.shade2) && !is.null(p$cint.upr)) {
            cints <- data.frame(upr=rinfo$scale * (p$cint.upr - fitted),
                                lwr=rinfo$scale * (p$cint.lwr - fitted))
            colnames(cints) <- c("upr", "lwr")
        }
    }
    if(is.null(pints) && is.null(cints))
        level <- 0

    resids <- rinfo$scale * rinfo$resids

    if((which %in% W7VLOG:W9LOGLOG))
        check.that.most.are.positive(
            versus1, "fitted", sprintf("which=%d", which), "nonpositive")

    # TODO following is redundant after above check?
    # abs(resids) must be nonnegative to take their log
    if(which %in% W7VLOG:W9LOGLOG)
        check.that.most.are.positive(
            abs(resids), "abs(residuals)", sprintf("which=%d", which), "zero")

    trans.versus <- trans.versus(versus1[iresids], which)
    trans.resids <- trans.resids(resids[iresids], which)
    x <- if(nversus == V2INDEX) 1:length(trans.versus) else trans.versus
    jitter <- as.numeric(check.numeric.scalar(jitter, logical.ok=TRUE))
    stopifnot(jitter >= 0, jitter <= 10) # 10 is arbitrary
    jittered.x <- x
    jittered.trans.resids <- trans.resids
    if(jitter > 0) {
        # we use amount=0 (same as S) which seems to work better in this context
        jittered.x            <- jitter(x,            factor=jitter, amount=0)
        jittered.trans.resids <- jitter(trans.resids, factor=jitter, amount=0)
    }
    derived.xlab <- derive.xlab(dot("xlab", DEF=NULL, ...),
                                which, colname.versus1, nversus)
    derived.ylab <- derive.ylab(dot("ylab", DEF=NULL, ...), which, rinfo$name)
    main <- derive.main(main=dot("main", DEF=NULL, ...),
                        derived.xlab, derived.ylab, level)

    # allow col.response as an argname for compat with old plotmo
    pt.col <- dot("col.response col.resp", DEF=1, ...)
    pt.col <- dot("pt.col col.points col.point col.residuals col.resid col",
                     EX=c(0,1,1,1,1,1), DEF=pt.col, NEW=1, ...)
    # recycle
    pt.col <- repl(pt.col, length(resids))

    pt.cex <- dot("response.cex cex.response", DEF=1, ...)
    pt.cex <- dot("pt.cex cex.points cex.point cex.residuals cex",
                     EX=c(0,1,1,1,1), DEF=pt.cex, NEW=1, ...)
    pt.cex <- pt.cex * pt.cex(length(x), npoints)
    pt.cex <- repl(pt.cex, length(resids))

    pt.pch <- dot("response.pch pch.response", DEF=20, ...)
    pt.pch <- dot("pt.pch pch.points pch.point pch.residuals pch",
                     EX=c(0,1,1,1,1), DEF=pt.pch, NEW=1, ...)
    pt.pch <- repl(pt.pch, length(resids))

    ylim <- get.resids.ylim(ylim=dot("ylim", ...), force.auto.resids.xylim,
                            object, fitted, trans.resids, which,
                            info, standardize, id.indices, center,
                            pints, cints, rinfo$scale, nversus)

    xlim <- get.resids.xlim(xlim=dot("xlim", ...), force.auto.resids.xylim,
                            which, x, trans.versus, ylim, nversus, id.indices)

    call.plot(graphics::plot.default, PREFIX="pt.",
         force.x    = x,
         force.y    = jittered.trans.resids,
         force.main = main,
         force.xlab = derived.xlab,
         force.ylab = derived.ylab,
         force.xlim = xlim,
         force.ylim = ylim,
         force.col  = NA, # no points will actually be plotted at this stage
         ...)
    if(is.specified(grid.col))
        grid(col=grid.col, lty=1)
    else if(which != W9LOGLOG)
        abline(h=0, lty=1, col="lightgray") # axis

    if(level && nversus != V4LEVER) {
        if(is.specified(level.shade))
            draw.pint.resids(pints=pints, x=versus1,
                             shade=level.shade, nversus=nversus, ...)
        if(is.specified(level.shade2))
            draw.pint.resids(pints=cints, x=versus1,
                             shade=level.shade2, nversus=nversus, ...)
    }
    if(nversus == V4LEVER) {
        # vertical line at mean leverage
        mean <- mean(x, na.rm=TRUE)
        abline(v=mean, col="gray")
        # add label "mean"
        if(which == W3RESID) { # not for others otherwise put text over the points
            usr <- par("usr") # xmin, xmax, ymin, ymax
            text(mean,
                 if(info) usr[3] + .1 * (usr[4] - usr[3]) # beyond density plot
                 else     usr[3] + .02 * (usr[4] - usr[3]),
                 "mean",
                 adj=c(0, -.2), cex=.8, srt=90)
        }
        if(standardize && inherits(object, "lm"))
            draw.cook.levels(object, ...)
    }
    call.plot(graphics::points, PREFIX="pt.",
              force.x   = jittered.x,
              force.y   = jittered.trans.resids,
              force.col = pt.col[iresids],
              force.cex = pt.cex[iresids],
              force.pch = pt.pch[iresids],
              ...)
    box()

    # plot points with unity leverage as stars
    draw.bad.leverage.as.star(jittered.x, rinfo, iresids, pt.cex, smooth.col)

    coef.rlm <- NULL
    if(info && nversus != V4LEVER && (which == W5ABS || which == W9LOGLOG))
        coef.rlm <- draw.rlm.line(which, versus1, resids, nversus, ...)

    if(which != W9LOGLOG)
        draw.smooth(x, trans.resids, rinfo$scale[iresids], smooth.col, ...)

    col.cv <- dot("col.cv", ...)
    oof.meanfit.was.plotted <- FALSE
    if(level && !is.null(object$cv.oof.fit.tab) && is.specified(col.cv))
    {
        draw.oof.meanfit(object$cv.oof.fit.tab, fitted, versus1, rinfo,
                         which, col.cv, nversus)
        oof.meanfit.was.plotted <- TRUE
    }
    # TODO implement id.indices for nversus=V2INDEX
    if(!is.null(id.indices) && nversus != V2INDEX) {
        # TODO as.numeric is needed if versus1 is a factor
        # is.na test needed for which=7 (if some are negative?)
        x1 <- as.numeric(trans.versus(versus1, which)[id.indices])
        if(!anyNA(x1))
            plotrix::thigmophobe.labels(x=x1,
                # TODO labels should take into account jitter
                y=trans.resids(resids, which)[id.indices],
                labels=rinfo$labs[id.indices],
                offset=.33, xpd=NA,
                font=dot("label.font", DEF=1, ...)[1],
                cex=.8 * dot("label.cex", DEF=1, ...)[1],
                col=dot("label.col",
                        DEF=if(is.specified(smooth.col)) smooth.col else 2, ...)[1])
    }
    if(info)
        draw.resids.info(which, info, versus1, resids, nversus, rsq, coef.rlm, ...)
    else
        possible.plotres.legend(which=which,
            level=level, smooth.col=smooth.col,
            oof.meanfit.was.plotted=oof.meanfit.was.plotted, ...)
}
get.plotres.data <- function(object, object.name, which, standardize, delever,
                             level, versus, id.n, labels.id,
                             trace, npoints, type, nresponse, ...)
{
    # the dot argument FORCEPREDICT is to check compat with old plot.earth
    meta <- plotmo_meta(object, type, nresponse, trace,
                        avoid.predict=!dot("FORCEPREDICT", DEF=FALSE, ...), ...)
      nresponse <- meta$nresponse  # column index
      resp.name <- meta$resp.name  # used only in automatic caption, may be NULL
      type      <- meta$type       # always a string (converted from NULL if necessary)
      residtype <- meta$residtype  # ditto

    rsq <- try(plotmo_rsq1(object=object, newdata=NULL,
                           trace=if(trace == 1) -1 else trace, meta=meta, ...),
                      silent=trace < 2)
    if(is.try.err(rsq))
        rsq <- NA

    # get the residuals and fitted info
    rinfo <- plotmo_rinfo(object=object, type=type, residtype=residtype,
                nresponse=nresponse,
                standardize=standardize, delever=delever, trace=trace,
                leverage.msg=
                    if(any(which %in% c(W3RESID,W5ABS:W9LOGLOG)))
                        "plotted as a star"
                    else
                        "ignored",
                expected.levs=meta$resp.levs, labels.id=labels.id, ...)

    fitted <- rinfo$fitted # n x 1 numeric matrix
    rinfo$fitted <- NA     # prevent accidental use of rinfo$fitted later

    stopifnot(NCOL(fitted) == 1)
    stopifnot(length(dim(fitted)) == 2)
    colnames(fitted) <- "Fitted" # colname will be used in labels in plots

    # get the values we will plot against (by default the fitted values)
    vinfo <- get.versus.info(which, versus, object, fitted, nresponse, trace)
    stopifnot(nrow(fitted) == length(rinfo$resids))

    ncases <- length(rinfo$resids)
    id.n <- get.id.n(id.n, ncases) # convert special values of id.n

    # convert special values of npoints
    check.integer.scalar(npoints, min=-1, null.ok=TRUE, logical.ok=TRUE)
    npoints.was.neg <- FALSE
    if(is.null(npoints))
        npoints <- 0
    else if(is.logical(npoints))
        npoints <- if(npoints) ncases else 0
    else if(npoints == -1) {
        npoints.was.neg <- TRUE
        npoints <- ncases
    } else if(npoints > ncases)
        npoints <- ncases

    # Use a maximum of NMAX residuals (unless npoints is bigger or negative).
    # Allows plotres to be fast even on models with millions of cases.
    NMAX <- 1e4
    nmax <- max(NMAX, npoints)
    if(!npoints.was.neg && nrow(fitted) > nmax) {
        if(trace >= 1)
            printf("using %g of %g residuals%s\n",
                nmax, ncases,
                if(id.n > 0)
                    ", forcing id.n=0 because of that (implementation restriction)"
                else
                    "")
        # see comment in plotres for use of V4LEVER here
        isubset <- get.isubset(rinfo$resids, nmax, id.n,
                     use.all=(vinfo$nversus == V4LEVER), rinfo$scale)
        fitted           <- fitted      [isubset, , drop=FALSE]
        rinfo$resids     <- rinfo$resids[isubset, , drop=FALSE]
        rinfo$scale      <- rinfo$scale [isubset]
        vinfo$versus.mat <- vinfo$versus.mat  [isubset, , drop=FALSE]
        # Can no longer draw point labels because row numbers are different.
        # TODO Come up with a solution so it doesn't have to be that way.
        id.n <- 0
    }
    list(nresponse  = nresponse, # col index in the response (converted from NA if necessary)
         resp.name  = resp.name, # used only in automatic caption, may be NULL
         type       = type,      # always a string (converted from NULL if necessary)
         rinfo      = rinfo,     # resids, scale, name, etc.
         vinfo      = vinfo,     # versus.mat, icolumns, nversus, etc.
         fitted     = fitted,    # n x 1 numeric matrix, colname is "Fitted"
         id.n       = id.n,      # forced to zero if row indexing changed
         npoints    = npoints,   # special values have been converted
         rsq        = rsq)
}
get.id.n <- function(id.n, ncases) # convert special values of id.n
{
    check.integer.scalar(id.n, null.ok=TRUE, logical.ok=TRUE)
    if(is.null(id.n))
        id.n <- 0
    else if(is.logical(id.n)) {
        id.n <- if(id.n) ncases else 0
    } else if(id.n == -1)
        id.n <- ncases
    else if(abs(id.n) > ncases)
        id.n <- ncases
    id.n
}
get.versus.info <- function(which, versus, object, fitted, nresponse, trace=0)
{
    versus.mat <- fitted
    icolumns   <- 1
    trim.which <- FALSE
    got.versus <- FALSE
    nversus    <- versus
    if(is.numeric(versus)) {
        got.versus <- TRUE
        trim.which <- TRUE
        if(length(versus) != 1)
            stop0(
"illegal 'versus' argument (length of 'versus' must be 1 when 'versus' is numeric)")
        if(floor(versus) != versus)
            versus.err()
        if(versus == V1FITTED)
            trim.which <- FALSE
        else if(versus == V2INDEX)
            NULL
        else if(versus == V3RESPONSE) {
            versus.mat <- plotmo_y(object, nresponse, trace,
                            expected.len=NROW(fitted), object$levels)$y
            colnames(versus.mat) <- "Response"
        } else if(versus == V4LEVER) {
            # TODO handle constant leverages for factors in the same way as plot.lm
            versus.mat <-
                matrix(hatvalues1(object, sprintf("versus=%g", V4LEVER)), ncol=1)
            colnames(versus.mat) <- "Leverage"
        } else
            versus.err()
    }
    else if(!is.character(versus))
        versus.err()
    else if(length(versus) == 1 && nchar(versus) >= 2 &&
            (substr(versus, 1, 2) == "b:" || substr(versus, 1, 2) == "B:")) {
        # use the basis matrix
        got.versus <- TRUE
        trim.which <- TRUE
        nversus <- 0
        plotmo_bx <- plotmo_bx(object, trace,
                               versus=substring(versus, 3)) # substring drops "bx:"
            versus.mat <- plotmo_bx$bx
            icolumns   <- plotmo_bx$icolumns
    }
    if(!got.versus) { # user specified x variables
        trim.which <- TRUE
        prefix <- substr(versus, 1, 1)
        nversus <- 0
        # following are needed if versus is a vector
        if(any(prefix == "*"))
            stop0("\"*\" is not allowed in this context in the 'versus' argument\n",
                  "      Your 'versus' argument is ", quote.with.c(versus))
        versus.mat <- plotmo_x(object, trace)
        versus.mat <- as.matrix(versus.mat)
        colnames(versus.mat) <- gen.colnames(versus.mat, "x", "x", trace)
        icolumns <- check.index(versus, "versus", seq_len(NCOL(versus.mat)),
                                colnames=colnames(versus.mat))
    }
    if(trim.which) {
        # remove all entries from which except standard resid and abs resid plots
        org.which <- which
        which <- which[which %in% c(W3RESID,W5ABS)]
        if(length(which) == 0)
            warnf(
                "which=%s is now empty because plots were removed because versus=%s",
                paste.c(org.which), paste(versus))
    }
    list(which      = which,      # which after possibly removing some plots
         versus.mat = versus.mat, # either fitted, response, x, or bx
         icolumns   = icolumns,   # desired column indices in versus.mat
         nversus    = nversus)    # versus as a number, 0 if versus is character
}
get.resids.xlim <- function(xlim, force.auto.resids.xylim,
                            which, x, trans.versus, ylim, nversus, id.indices)
{
    if(force.auto.resids.xylim || !is.specified(xlim)) { # auto xlim?
        if(which == W9LOGLOG) {
            # don't show lower 5% of points
            quant <- quantile(trans.versus, prob=c(.05, 1), na.rm=TRUE)
            min <- quant[1]
            max <- quant[2]
            # extra left margin so slope of linear fit not flattened
            if(min > .2 * ylim[1])
                min <- .2 * ylim[1]
            xlim <- c(min, max)
        } else if(nversus == V4LEVER) # room for labels on high leverage points
            xlim <- c(0, 1.1 * max(x, na.rm=TRUE))
        else
            xlim <- range1(x, na.rm=TRUE)
        range <- xlim[2] - xlim[1]
        if(is.specified(id.indices)) # space for point labels
            xlim <- c(xlim[1] - .04 * range, xlim[2] + .04 * range)
    }
    stopifnot(is.numeric(xlim), length(xlim) == 2)
    fix.lim(xlim)
}
get.resids.ylim <- function(ylim, force.auto.resids.xylim,
                            object, fitted, resids, which,
                            info, standardize, id.indices, center,
                            pints, cints, scale, nversus)
{
    if(force.auto.resids.xylim || !is.specified(ylim)) { # auto xlim?
        if(!is.null(pints)) {
            min <- min(resids, pints$lwr, na.rm=TRUE)
            max <- max(resids, pints$upr, na.rm=TRUE)
        } else if(!is.null(cints)) {
            min <- min(resids, cints$lwr, na.rm=TRUE)
            max <- max(resids, cints$upr, na.rm=TRUE)
        } else {
            min <- min(resids, na.rm=TRUE)
            max <- max(resids, na.rm=TRUE)
        }
        maxa <- mina <- 0 # adjustments to max and min
        if(which %in% (W5ABS:W8CUBE))
            min <- 0
        else if(which == W3RESID && center) {
            # want symmetric ylim so can more easily see asymmetry
            if(abs(min) > abs(max))
                max <- -min
            else if(abs(max) > abs(min))
                min <- -max
        } else if(which == W9LOGLOG)
            maxa <- .5  # more space on top, looks better
        range <- abs(max - min)
        if(is.specified(id.indices)) { # space for point labels
            # TODO only do this if point labels are near the edges
            mina <- max(mina, .03 * range)
            maxa <- max(maxa, .03 * range)
        }
        if(nversus == V4LEVER && standardize && inherits(object, "lm")) {
            maxa <- max(maxa, maxa + .2 * range) # space for cook distance legend
            mina <- max(mina, mina + .1 * range) # space for "mean" label
        }
        if(info) { # space for extra text (on top) and density plot (in the bottom)
            maxa <- maxa + max * if(!is.null(id.indices)) .2 else .1
            mina <- mina + max * if(!is.null(id.indices)) .2 else .1
        }
        ylim <- c(min-mina, max+maxa)
    }
    fix.lim(ylim)
}
draw.pint.resids <- function(pints, x, shade, nversus, ...)
{
    if(!is.null(pints)) {
        # abscissa must be ordered for polygon to work
        order <- order(x)
        x <- x[order]
        pints <- pints[order,]
        x <- if(nversus == V2INDEX)
                c(1:length(x), length(x):1)
             else
                trans.versus(c(x, rev(x)), 0)
        call.plot(graphics::polygon, PREFIX="level.", drop.shade=1, drop.shade2=1,
            force.x    = x,
            force.y    = trans.resids(c(pints$lwr, rev(pints$upr)), 0),
            force.col  = shade,
            def.border = shade,
            def.lty    = 0,
            ...)
    }
}
# this should be used only for models with homoscedastic errors
draw.cook.levels <- function(object, ...)
{
    cook.levels <- dot("cook.levels", DEF=c(0.5, 1.0), ...)
    stopifnot(is.numeric(cook.levels), all(cook.levels > 0))
    col <- dot("cook.col", DEF="slategray4", ...)
    lty <- dot("cook.lty", DEF=1,          ...)
    lwd <- dot("cook.lwd", DEF=1,          ...)
    # based on code in stats::plot.lm.R
    leverage <- hatvalues1(object, "'standardize'")
    p <- length(coef(object))
    leverage.range <- range(leverage, na.rm=TRUE) # though should never have NA
    x <- seq.int(0, 1, length.out=101)
    for(cook.level in cook.levels) {
        cl <- sqrt(cook.level * p *(1 - x) / x)
        lines(x,  cl, col=col, lty=lty, lwd=lwd)
        lines(x, -cl, col=col, lty=lty, lwd=lwd)
    }
    # we don't use bottomleft like plot.lm because we may plot the density there
    usr <- par("usr") # xmin, xmax, ymin, ymax
    legend(usr[1]-.7 * strwidth("X"), # jam it into the corner
           usr[4]+.5 * strheight("X"),
           legend="Cook's distance", col=col, lty=lty, lwd=lwd,
           box.col="white", bg="white", x.intersp=.2, seg.len=1.5)
    xmax <- min(0.99, usr[2])
    ymult <- sqrt(p * (1 - xmax) / xmax)
    axis(4, at=c(-sqrt(rev(cook.levels)) * ymult, sqrt(cook.levels)*ymult),
         labels=paste(c(rev(cook.levels), cook.levels)),
         mgp=c(.25,.15,0), las=2, tck=0,
         cex.axis=.7, col.axis=col,
         font=2) # makes the gray labels a bit more legible
}
# Plot points with unity leverage as stars.  We plot them on
# the axis, which is arguably incorrect but still useful.
# TODO add a test for this to the test suite

draw.bad.leverage.as.star <- function(x, rinfo, iresids, pt.cex, smooth.col)
{
    which <- which(is.na(rinfo$scale[iresids]))
    if(length(which) > 0) {
        points(x[which], 0, col=1, cex=pt.cex[iresids], pch=8) # pch 8 is a star
        # add label if possible (not poss if not all points plotted, see npoints)
        if(length(iresids) == length(rinfo$scale)) {
            label <- which(is.na(rinfo$scale))
            text.on.white(x=x[which], y=0, label=label,
                          col=if(is.specified(smooth.col)) smooth.col else 2,
                          cex=.8, adj=-.5, xpd=NA)
        }
    }
}
draw.smooth <- function(x, resids, scale, smooth.col, ...)
{
    if(!is.specified(smooth.col))
        return(NULL)

    # na.rm is needed if we take logs of nonpos, see check.that.most.are.positive.
    # That's why we calculate delta explicitly instead of using lowess default.
    delta <- .01 * diff(range1(x, na.rm=TRUE))

    # Replace points with NA scale with 0 (else lowess stops at the NA).
    # Zero is appropriate because the points are 0 resids with leverage 1.
    resids[which(is.na(scale))] <- 0

    # we use lowess rather than loess because loess tends to give warnings
    smooth.f    <- dot("smooth.f loess.f", DEF=2/3, NEW=1, ...)
    smooth.iter <- dot("smooth.iter",      DEF=3,          ...)
    check.numeric.scalar(smooth.f)
    stopifnot(smooth.f > .01, smooth.f < 1)
    smooth <- lowess(x, resids, f=smooth.f, iter=smooth.iter, delta=delta)
    call.plot(graphics::lines.default, PREFIX="smooth.", drop.f=1,
          force.x   = smooth$x,
          force.y   = smooth$y,
          force.col = smooth.col,
          force.lwd = dot("smooth.lwd lwd.smooth lwd.loess",
                          EX=c(0,1,1), DEF=1, NEW=1, ...),
          force.lty = dot("smooth.lty lty.smooth",
                          EX=c(0,1), DEF=1, NEW=1, ...),
          ...)
}
derive.xlab <- function(xlab, which, colname.versus1, nversus)
{
    if(is.specified(xlab)) {
        stopifnot.string(xlab, allow.empty=TRUE)
        if(!nzchar(xlab))
            return("")
    }
    if(!is.specified(xlab))
        xlab <- colname.versus1
    stopifnot.string(xlab)
    if(which %in% (W7VLOG:W9LOGLOG))
        xlab <- sprintf("Log %s", xlab)
    if(nversus == V2INDEX)
        xlab <- sprintf("%s index", xlab)
    xlab
}
derive.ylab <- function(ylab, which, rinfo.name)
{
    if(is.specified(ylab)) {
        stopifnot.string(ylab, allow.empty=TRUE)
        if(!nzchar(ylab))
            return("")
    }
    if(!is.specified(ylab))
        ylab <- sprintf("%ss", rinfo.name)
    if(which == W5ABS)
        ylab <- sprintf("Abs %s", ylab)
    else if(which == W6SQRT)
        ylab <- sprintf("Sqrt Abs %s", ylab)
    else if(which == W7VLOG)
        ylab <- sprintf("Abs %s", ylab)
    else if(which == W8CUBE)
        ylab <- sprintf("Cube Root Squared %s", ylab)
    else if(which == W9LOGLOG)
        ylab <- sprintf("Log Abs %s", ylab)
    ylab
}
derive.main <- function(main, xlab, ylab, level) # title of plot
{
    # TODO should really use strwidth for newline calculation
    # The "Fitted" helps with limitations of the formula below
    newline <- xlab != "Fitted" && xlab != "Fitted index" && xlab != "Response" &&
               nchar(ylab) + nchar(xlab) > 15
    if(xlab == "Leverage" && ylab == "Residuals") # special case, mainly for which=1 with lm
        newline <- FALSE
    else if(grepl("Standardized", ylab[1]) || grepl("Delevered", ylab[1]))
        newline <- TRUE

    if(!is.specified(main)) # generate a main only if user didn't specify main
        main <- sprintf("%s vs%s%s", ylab, if(newline) "\n" else " ", xlab)
    if(xlab != "Leverage" && level && !newline) # two newlines is too many
        main <- sprintf("%s\n%g%% level shaded", main, 100*(level))

    main
}
# plot resids of oof meanfit with col.cv (default lightblue)

draw.oof.meanfit <- function(cv.oof.fit.tab, fitted, versus1,
                             rinfo, which, col.cv, nversus)
{
    # mean of each row of oof.fit.tab
    meanfit <- apply(cv.oof.fit.tab,  1, mean)
    meanfit <- rinfo$scale * (meanfit - fitted)
    order <- order(versus1)
    trans.versus1 <- trans.versus(versus1[order], which)
    x <- if(nversus == V2INDEX) 1:length(trans.versus1) else trans.versus1
    lines(x, trans.resids(meanfit[order], which), col=col.cv)
}
draw.density.along.the.bottom <- function(x, den.col=NULL, scale=NULL, ...)
{
    if(is.null(den.col))
        den.col <- dot("density.col", DEF="gray57", EX=0, ...)
    den <- density(x, adjust=dot("density.adjust", DEF=.5, EX=0, ...), na.rm=TRUE)
    usr <- par("usr") # xmin, xmax, ymin, ymax
    if(is.null(scale))
        scale <- .08 / (max(den$y) - min(den$y))
    den$y <- usr[3] + den$y * scale * (usr[4] - usr[3])
    call.plot(graphics::lines.default, PREFIX="density.", drop.adjust=1,
              force.x=den, force.y=NULL, def.col=den.col, ...)
}
draw.rlm.line <- function(which, versus1, resids, nversus, ...)
{
    trans.resids <- trans.resids(resids, which)
    trans.versus <- trans.versus(versus1, which)
    x <- if(nversus == V2INDEX) 1:length(trans.versus) else trans.versus
    if(which == W9LOGLOG) {
        # ignore lower 10% of points (very small residuals i.e. very neg logs)
        quant <- quantile(trans.versus, prob=.1, na.rm=TRUE)
        ok <- which(x > quant)
        x <- x[ok]
        trans.resids <- trans.resids[ok]
    }
    # # regression on 10 bootstrap samples so we can see variance of versus1
    # for(i in 1:10) {
    #   j <- sample.int(length(x), replace=TRUE)
    #   trans.resids1 <- trans.resids[j]
    #   trimmed.trans.fit1 <- x[j]
    #   rlm <- MASS::rlm(trans.resids1~trimmed.trans.fit1,
    #                    method="MM", na.action="na.omit")
    #   if(draw)
    #       abline(coef(rlm), col="lightgray", lwd=.6)
    # }

    # robust linear regression of trans.resids on x
    # na.omit is needed if some versus1 (or resids) were nonpos so log(versus1) is NA
    rlm <- MASS::rlm(trans.resids~x, method="MM", na.action="na.omit")
    call.dots(abline,
        force.coef = coef(rlm),
        force.col  = "lightblue",
        force.lwd  = dot("smooth.lwd lwd.smooth lwd.loess",
                         EX=c(0,1,1), DEF=1, NEW=1, ...) + 1,
        ...)
    box() # abline overplots the box
    coef(rlm)
}
draw.resids.info <- function(which, info, versus1, resids, nversus, rsq, coef.rlm, ...)
{
    trans.versus <- trans.versus(versus1, which)
    x <- if(nversus == V2INDEX) 1:length(trans.versus) else trans.versus
    # TODO consider also drawing the density along the right side
    draw.density.along.the.bottom(x, ...)
    if(nversus != V4LEVER) {
        lm.text <- ""
        slope.text <- ""
        if(which == W5ABS || which == W9LOGLOG) { # added linear regression line?
            stopifnot(length(coef.rlm) == 2)
            slope.text <- sprintf("    slope %.2g", coef.rlm[2])
        }
        # exact=FALSE else get warning "Cannot compute exact p-value with ties"
        cor.abs <- cor.test(versus1, abs(resids), method="spearman", exact=FALSE)
        if(nversus == V3RESPONSE) {
            cor <- cor.test(versus1, resids, method="spearman", exact=FALSE)
            text <- sprintf("spearman abs  %.2f   resids %.2f\n%s",
                             cor.abs$estimate, cor$estimate, slope.text)
        } else if(which == W3RESID && nversus == V1FITTED)
            text <- sprintf("rsq  %.2f     spearman abs  %.2f",
                            rsq, cor.abs$estimate)
        else
            text <- sprintf("spearman abs  %.2f%s", cor.abs$estimate, slope.text)
        cex <- .9
        usr <- par("usr") # xmin, xmax, ymin, ymax
        text.on.white(x     = usr[1] + strwidth("x", font=1),
                      y     = usr[4] - cex * (strheight(text, font=1) +
                                  .5 * strheight("X", font=1)),
                      label = text,
                      cex   = cex, adj=c(0, 0), font=1, col=1, xpd=NA)
    }
}
my.log10 <- function(x) # log of very small values is silently set to NA
{
    i <- which(x < max(x) / 1e6)
    x[i] <- 1
    x <- log10(x)
    x[i] <- NA
    x
}
trans.versus <- function(versus, which)
{
    if(which %in% (W7VLOG:W9LOGLOG))
        my.log10(versus)
    else
        versus
}
trans.resids <- function(resid, which) # transform the residuals
{
    if(which == W5ABS)
        abs(resid)
    else if(which == W6SQRT)
        sqrt(abs(resid))
    else if(which == W7VLOG)
        abs(resid)
    else if(which == W8CUBE) {
        # do it in two steps so no problem with negative numbers
        resid <- resid^2
        resid^(1/3)
    } else if(which == W9LOGLOG)
        my.log10(abs(resid))
    else
        resid
}
# Get a subset of x.  Size of subset is nsubset.  Returns an index vector.
# The subset includes the 20 biggest absolute values in x.
# If you want more than the 20 biggest values, set nbiggest.

get.isubset <- function(x, nsubset, nbiggest=0, use.all=FALSE, scale=NULL)
{
    check.vec(x, "x passed to get.isubset", length(x))
    ix <- seq_len(length(x))
    if(nsubset > 0 && nsubset < length(x) && !use.all) { # TODO move this into caller
        # take a sample, but make sure it includes the 20 biggest absolute values
        nsubset <- min(nsubset, length(x))
        nbiggest <- min(max(20, nbiggest), nsubset)
        isorted <- order(abs(x), decreasing=TRUE) # expensive
        ikeep <- seq_len(nbiggest)
        if(nsubset > nbiggest)
            ikeep <- c(ikeep, seq(from=nbiggest + 1, to=length(x),
                                  length.out=nsubset - nbiggest))
        ix <- isorted[ikeep]
        # force in points with unity leverage
        if(!is.null(scale)) {
            which <- which(is.na(scale))
            if(length(which) > 0)
                ix <- sort.unique(c(which, ix))
        }
    }
    ix # index vector
}
# get the indices of the id.n biggest resids (requires a sort)

get.id.indices <- function(resids, id.n, hatvalues=NULL)
{
    # id.n has already been checked in plotres.data
    if(id.n == 0)
        return(NULL)
    if(id.n > 0) { # same as plot.lm
        i <- sort.list(abs(resids), decreasing=TRUE, na.last=NA)
        if(length(i) > id.n)
            i <- i[1:id.n]
    } else { # id.n is negative: most positive and most negative residuals
        id.n <- -id.n
        i <- sort.list(resids, decreasing=TRUE, na.last=NA)
        if(length(i) > id.n)
            i <- i[c(1:id.n, length(i) + 1 - (1:id.n))]
    }
    if(!is.null(hatvalues)) {
        # add the worst hatvalues i.e. the worst leverages
        i <- unique(c(i, order(hatvalues, decreasing=TRUE)[1:id.n]))
    }
    i
}
possible.plotres.legend <- function(which, level, smooth.col,
                                    oof.meanfit.was.plotted, ...)
{
    # add legend, else red and blue may confuse the user
    legend.pos <- dot("legend.pos", DEF=NULL, ...)
    if(level && oof.meanfit.was.plotted &&
           (is.null(legend.pos) || !all(is.na(legend.pos)))) {
        if(is.null(legend.pos)) { # auto?
            legend.x <- "bottomleft"
            legend.y <- NULL
        } else { # user specified legend position
            legend.x <- legend.pos[1]
            legend.y <- if(length(legend.pos) > 1) legend.pos[2] else NULL
        }
        legend.txt <- NULL
        legend.col <- NULL
        legend.lwd <- NULL
        legend.lty <- NULL
        if(which != W9LOGLOG && is.specified(smooth.col)) { # smooth plotted?
            legend.txt <- "smooth"
            legend.col <- smooth.col
            legend.lwd <- dot("smooth.lwd lwd.smooth lwd.loess",
                              EX=c(0,1,1), DEF=1, NEW=1, ...)
            legend.lty <- 1
        }
        if(oof.meanfit.was.plotted) {
            legend.txt <- c(legend.txt, "cross validated oof fit")
            legend.col <- c(legend.col, dot("col.cv", ...))
            legend.lwd <- c(legend.lwd, 1)
            legend.lty <- c(legend.lty, 1)
        }
        if(!is.null(legend.txt))
            call.dots(graphics::legend,
                DROP="*", KEEP="PREFIX",
                force.x      = legend.x,
                force.y      = legend.y,
                force.legend = legend.txt,
                def.col      = legend.col,
                def.lwd      = legend.lwd,
                def.lty      = legend.lty,
                def.bg       = "white",
                def.cex      = .8,
                ...)
    }
}
