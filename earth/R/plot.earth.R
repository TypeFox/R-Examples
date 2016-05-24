# plot.earth.R: plotting routines for the earth package

plot.earth <- function(x = stop("no 'x' argument"),
    which       = 1:4,
    info        = FALSE,
    versus      = 1,

    standardize = FALSE,
    delever     = FALSE,
    level       = 0,

    id.n        = 3,
    labels.id   = NULL,
    smooth.col  = 2,
    grid.col    = 0,
    jitter      = 0,

    do.par      = NULL,
    caption     = NULL,
    trace       = 0,

    npoints     = 3000,
    center      = TRUE,

    type        = NULL, # passed to predict
    nresponse   = NA,

    # following are passed to plotres via plotres's dots

    col.cv       = "lightblue",

    # following are passed to earth_plotmodsel via plotres's dots

    col.grsq            = 1,
    col.rsq             = 2,
    col.infold.rsq      = 0,
    col.mean.infold.rsq = 0,
    col.mean.oof.rsq    = "palevioletred",
    col.npreds          = if(is.null(object$cv.oof.rsq.tab)) 1 else 0,
    col.oof.labs        = 0,
    col.oof.rsq         = "mistyrose2",
    col.oof.vline       = col.mean.oof.rsq,
    col.pch.cv.rsq      = 0,
    col.pch.max.oof.rsq = 0,
    col.vline           = col.grsq,
    col.vseg            = 0,
    lty.grsq            = 1,
    lty.npreds          = 2,
    lty.rsq             = 5,
    lty.vline           = "12",

    legend.pos          = NULL,

    ...)
{
    object.name <- quote.deparse(substitute(x))
    object <- x
    remove(x) # prevent confusion with the x matrix
    check.classname(object, substitute(object), "earth")

    npoints  <- dot("nresiduals", DEF=npoints,  ...) # back compat
    col.rsq  <- dot("col.line",   DEF=col.rsq,  ...)

    plotmo::plotres(object=object, which=which, info=info, versus=versus,
        standardize=standardize, delever=delever, level=level,
        id.n=id.n, labels.id=labels.id, smooth.col=smooth.col,
        grid.col=grid.col, jitter=jitter,
        do.par=do.par, caption=caption, trace=trace,
        npoints=npoints, center=center,
        # type is passed in dots, it's not really a plot.earth arg
        nresponse=nresponse,
        object.name=object.name,
        # following are passed to plotres via plotres's dots
        col.cv=col.cv,
        # following are passed to earth_plotmodsel
        w1.col.grsq            = col.grsq,
        w1.col.rsq             = col.rsq,
        w1.col.infold.rsq      = col.infold.rsq,
        w1.col.mean.infold.rsq = col.mean.infold.rsq,
        w1.col.mean.oof.rsq    = col.mean.oof.rsq,
        w1.col.npreds          = col.npreds,
        w1.col.oof.labs        = col.oof.labs,
        w1.col.oof.rsq         = col.oof.rsq,
        w1.col.oof.vline       = col.oof.vline,
        w1.col.pch.cv.rsq      = col.pch.cv.rsq,
        w1.col.pch.max.oof.rsq = col.pch.max.oof.rsq,
        w1.col.vline           = col.vline,
        w1.col.vseg            = col.vseg,
        w1.lty.grsq            = lty.grsq,
        w1.lty.npreds          = lty.npreds,
        w1.lty.rsq             = lty.rsq,
        w1.lty.vline           = lty.vline,
        w1.legend.pos          = legend.pos,
        ...)
}
# TODO add nresponse to plot.earth.models

plot.earth.models <- function(
    x            = stop("no 'x' argument"),
    which        = c(1:2),
    caption      = "",
    jitter       = 0,
    col.grsq     = discrete.plot.cols(length(objects)),
    lty.grsq     = 1,
    col.rsq      = 0,
    lty.rsq      = 5,
    col.vline    = col.grsq,
    lty.vline    = "12",
    col.npreds   = 0,
    lty.npreds   = 2,
    legend.text  = NULL,
    do.par       = NULL,
    trace        = 0,
    ...)
{
    objects <- x
    remove(x) # prevent confusion with the x matrix
    if(!is.list(objects))       # note that is.list returns TRUE for a single object
        stop0("'x' is not an \"earth\" object or a list of \"earth\" objects")
    trace <- as.numeric(check.integer.scalar(trace, logical.ok=TRUE))
    # check for a common error, using plot.earth.models(mod1, mod2) instead
    # of plot.earth.models(list(mod1, mod2)) instead
    if(inherits(which, "earth"))
        stop0("use plot.earth.models(list(model1, model2)), ",
              "not plot.earth.models(model1, model2)")
    if(typeof(objects[[1]]) != "list") # if user specified just one object, convert to list
        objects <- list(objects)
    check.index(which, "which", 1:2)
    show <- to.logical(which, 4)
    if(length(which) == 0) {
        warning0("plot.earth.models: nothing to plot (the 'which' argument is empty)")
        return(invisible())
    }
    if(is.null(col.rsq))
        col.rsq <- if(is.null(col.grsq)) col.rsq else col.grsq
    if(is.null(col.npreds))
        col.npreds <- if(is.null(col.grsq)) col.rsq else col.grsq
    cum.col1 <- dot("cum.col col.cum pt.col col", ...)
    if(!is.specified(cum.col1))
        cum.col1 <- if(!is.specified(col.grsq)) col.rsq else col.grsq
    if(show[1] && col.grsq[1] == 0 && col.rsq[1] == 0)
        stop0("both col.grsq[1] and col.rsq[1] are zero")
    if(show[2] && !is.specified(cum.col1))
        stop0("cum.col is NULL, and unable to use col.grsq or col.rsq instead")
    nmodels <- length(objects)
    col.grsq   <- repl(col.grsq, nmodels)
    lty.grsq   <- repl(lty.grsq, nmodels)
    col.rsq    <- repl(col.rsq, nmodels)
    lty.rsq    <- repl(lty.rsq, nmodels)
    col.npreds <- repl(col.npreds, nmodels)
    lty.npreds <- repl(lty.npreds, nmodels)
    cum.col1   <- repl(cum.col1, nmodels)
    col.vline  <- repl(col.vline, nmodels)
    lty.vline  <- repl(lty.vline, nmodels)
    do.par <- check.do.par(do.par, length(which)) # do.par is 0, 1, or 2
    # prepare caption --- we need it now for do.par() but
    # can only display it later after at least one plot
    if(is.null(caption))
        caption <- ""
    main <- dot("main", DEF="Model Comparison", ...)
    if(do.par) {
        oldpar <- par(no.readonly=TRUE)
        do.par(nfigs=length(which), caption=caption, main1=main,
               xlab1=NULL, ylab1=NULL, trace=trace,
               def.font.main=1, ...) # for compat with lm.plot
        if(do.par == 1)
            on.exit(par(oldpar), add=TRUE)
    } else { # do.par=FALSE
        oldpar <- do.par.dots(..., trace=trace)
        if(length(oldpar))
            on.exit(do.call(par, oldpar), add=TRUE)
    }
    max.npreds <- 1
    max.nterms <- 1
    ylim <- dot("ylim", DEF=c(0,1), ...)
    for(imodel in seq_along(objects)) {
        object <- objects[[imodel]]
        check.classname(object, objects[[imodel]], "earth")
        ylim <- range(ylim,
                      get.model.selection.ylim(object, ylim,
                                               col.grsq[imodel], col.rsq[imodel]))
        max.npreds <- max(max.npreds,
                          get.nused.preds.per.subset(object$dirs, object$prune.terms))
        max.nterms <- max(max.nterms, length(object$rss.per.subset))
    }
    legend.col <- dot("legend.col col.legend", EX=c(0,1), DEF=1, NEW=1, ...)
    if(show[1]) {
        if(is.null(object$residuals)) # probably a model from object$cv.list
            stop0("earth object has no $residuals field.\n",
                  "       Use keepxy=TRUE in the call to earth.")
        for(imodel in seq_along(objects))
            earth_plotmodsel(
                x                   = objects[[imodel]],
                col.rsq             = col.rsq[imodel],
                col.grsq            = col.grsq[imodel],
                col.infold.rsq      = 0,
                col.mean.infold.rsq = 0,
                col.mean.oof.rsq    = 0,
                col.npreds          = col.npreds[imodel],
                col.oof.labs        = 0,
                col.oof.rsq         = 0,
                col.oof.vline       = 0,
                col.pch.cv.rsq      = 0,
                col.pch.max.oof.rsq = 0,
                col.vline           = col.vline[imodel],
                col.vseg            = col.grsq[imodel],
                lty.grsq            = lty.grsq[imodel],
                lty.npreds          = lty.npreds[imodel],
                lty.rsq             = lty.rsq,
                lty.vline           = lty.vline[imodel],

                legend.pos          = NA, # we plot our own legend
                add                 = (imodel > 1),
                max.nterms          = max.nterms,
                max.npreds          = max.npreds,
                jitter              = if(imodel>1) jitter else 0,

                # dots args
                main                = if(imodel > 1) "" else main,
                ylim                = ylim)

        if(is.specified(legend.col) && length(objects) > 1 && !show[2])
            draw.earth.models.legend(objects, min.width=.4,
                legend.text, legend.col,
                col.rsq, lty.rsq, col.grsq, lty.grsq, ...)
    }
    if(show[2]) {
        multiple.responses <- FALSE
        xlim <- c(0,0)
        for(object in objects) {
            if(is.null(object$residuals)) # probably a model from object$cv.list
                stop0("earth object has no $residuals field.\n",
                      "       Use keepxy=TRUE in the call to earth.")
            if(NCOL(object$residuals) > 1) {
                multiple.responses <- TRUE
                xlim <- range(xlim, abs(object$residuals[,1]), na.rm=TRUE)
            } else
                xlim <- range(xlim, abs(object$residuals), na.rm=TRUE)
        }
        for(imodel in seq_along(objects)) {
            object <- objects[[imodel]]
            attr(object, ".Environment") <- get.model.env(object)
            rinfo <- plotmo::plotmo_rinfo(object, type="earth", residtype="earth",
                                          nresponse=NULL,
                                          trace=trace, leverage.msg="ignored")
            plotmo::plotmo_cum(
                rinfo    = rinfo,
                info     = FALSE,
                nfigs    = 1,
                add      = (imodel > 1),
                cum.col1 =
                    if(length(cum.col1) > 1)                cum.col1[imodel]
                    else if(is.specified(col.grsq[imodel])) col.grsq[imodel]
                    else                                    col.rsq[imodel],
                grid.col = 0,
                jitter   = if(imodel == 1) 0 else jitter,
                cum.grid = "none",
                # dots args
                xlim     = xlim,
                main     =
                    if(imodel > 1)              ""
                    else if(multiple.responses) "Cumul Distrib (response 1)"
                    else                        "Cumulative Distribution")
        }
        if(is.specified(legend.col) && length(objects) > 1)
            draw.earth.models.legend(objects, min.width=.5,
                legend.text, legend.col,
                col.rsq, lty.rsq, col.grsq, lty.grsq, ...)
    }
    draw.caption(caption, ...)
    invisible()
}
# Return a vector of n clearly distinguishable colors.
# The first three are also distinguishable on (my) monochrome printer.

discrete.plot.cols <- function(ncolors=5)
{
    cols <- c(1, "brown", "gray60", "lightblue", "pink", "green")
    if(ncolors > length(cols))   # won't really be distinguishable
        cols <- c(cols, heat.colors(ncolors - length(cols)))
    cols[seq_len(ncolors)]
}
draw.earth.models.legend <- function(
    objects,
    min.width,
    legend.text,
    legend.col,
    col.rsq,
    lty.rsq,
    col.grsq,
    lty.grsq,
    ...)
{
    lty <- NULL
    col <- NULL
    if(is.null(legend.text)) {
        if(is.null(names(objects))) {
            args <- get.arg.strings(objects, maxchars=20)
            legend.text <- character(length=length(objects))
            for(imodel in seq_along(objects))
                legend.text[imodel] <- paste(imodel, args[[imodel]])
        } else
            legend.text <- names(objects)
    } else
        legend.text <- repl(legend.text, length(objects))
     if(col.rsq[1] != 0) {       # RSq plotted?
        col <- c(col, col.rsq)
        lty <- c(lty, repl(lty.rsq, length(col)))
        if(col.grsq[1] != 0)
            legend1 <- paste("RSq", legend.text)
    }
    if(col.grsq[1] != 0) {      # GRSq plotted?
        col <- c(col, col.grsq)
        lty <- c(lty, repl(lty.grsq, length(col)))
        if(col.rsq[1] != 0)
            legend.text <- c(legend1, paste("GRSq", legend.text))
    }
    legend.pos <- dot("legend.pos", DEF=NULL, ...)
    if(is.null(legend.pos)) { # auto?
        legend.x <- "bottomright"
        legend.y <- NULL
    } else { # user specified legend position
        legend.x <- legend.pos[1]
        legend.y <- if(length(legend.pos) > 1) legend.pos[2] else NULL
    }
    legend.cex <- get.earth.legend.cex(legend.text, min.width=min.width, ...)
    elegend(x=legend.x, y=legend.y, bg="white",
            legend=legend.text, col=col, lty=lty, cex=legend.cex,
            # y offset allows vertical lines to be visible below legend
            inset=c(.02, .04))
}
# called by plotres for which=1, and called by plot.earth.models
earth_plotmodsel <- function(
    x,
    col.rsq             = 2,
    col.grsq            = 1,
    col.infold.rsq      = 0,
    col.mean.infold.rsq = 0,
    col.mean.oof.rsq    = "palevioletred",
    col.npreds          = NULL,
    col.oof.labs        = 0,
    col.oof.rsq         = "mistyrose2",
    col.oof.vline       = col.mean.oof.rsq,
    col.pch.cv.rsq      = 0,
    col.pch.max.oof.rsq = 0,
    col.vline           = col.grsq,
    col.vseg            = 0,
    lty.grsq            = 1,
    lty.npreds          = 2,
    lty.rsq             = 5,
    lty.vline           = "12",

    legend.pos          = NULL,
    add        = FALSE,
    jitter     = 0,
    max.nterms = length(object$rss.per.subset),
    max.npreds = max(1,
                     get.nused.preds.per.subset(object$dirs, object$prune.terms)),

    ...)
{
    possibly.issue.cv.warning <- function()
    {
        if((!identical(col.mean.oof.rsq, "palevioletred") && !identical(col.mean.oof.rsq, 0)) ||
           (!identical(col.oof.rsq, "mistyrose2")         && !identical(col.oof.rsq, 0))      ||
           !identical(col.oof.labs, 0)        ||
           !identical(col.pch.max.oof.rsq, 0) ||
           !identical(col.pch.cv.rsq, 0)      ||
           !identical(col.mean.infold.rsq, 0) ||
           !identical(col.infold.rsq, 0)) {
            # user specifed a cross-validation argument, check that data is available
            if(is.null(object$cv.list))
                warning0("no cross-validation data because nfold not used in original call to earth")
            else if(is.null(object$cv.oof.rsq.tab))
                warning0("cannot plot cross-validation data because ",
                         "the earth model was not built with keepxy=TRUE")
        }
    }
    scale1 <- function(x, Min, Max)
    {
        return((x-Min)/(Max-Min))
    }
    left.axis <- function()
    {
        pretty <- pretty(c(ylim[1], ylim[2]))
        axis(side=2, at=scale1(pretty, ylim[1], ylim[2]), labels=pretty, srt=90)
        text <- ""
        if(is.specified(col.grsq))
            text <- "GRSq"
        if(is.specified(col.rsq) ||
                is.specified(col.oof.rsq) || is.specified(col.mean.oof.rsq) ||
                is.specified(col.infold.rsq) || is.specified(col.mean.infold.rsq))
            text <- paste0(text, "   RSq")
        # TODO mtext needs cex=par("cex"), not sure why
        # the line setting depends on the axis margin lines (want
        # compact axes if do.par set, but not compact if not set)
        mtext(text, side=2, cex=par("cex"),
              line=if(par("mgp")[1] < 1.8) 1.6 else 2)
    }
    right.axis <- function()
    {
        if(max.npreds <= 5) # try to get rid of fractions in the label
            pretty <- pretty(c(0, max.npreds), n=max.npreds)
        else
            pretty <- pretty(c(0, max.npreds))
        axis(side=4, at=scale1(pretty, 0, max.npreds), labels=pretty, srt=90)
        mtext("Number of used predictors", side=4, cex=par("cex"),
              line=if(par("mgp")[1] < 1.8) 1.4 else 1.8)
    }
    draw.selection.grid <- function() # plot the grid
    {
        if(!is.specified(grid.col))
            return()
        col <- grid.col[1]
        abline(v=0:par("usr")[2], col=col) # vertical grid
        if((ylim[2] - ylim[1]) > .5) # coarse horizontal grid?
            for(v in seq(-1, 1, by=.05))
                abline(h=scale1(v, ylim[1], ylim[2]), col=col, lwd=1)
        else { # fine horizontal grid
            for(v in seq(-1, 1, by=.01))
                abline(h=scale1(v, ylim[1], ylim[2]), col=col, lwd=.6)
            for(v in seq(-1, 1, by=.05))
                abline(h=scale1(v, ylim[1], ylim[2]), col=col, lwd=1.2)
        }
    }
    draw.infold.rsqs <- function() # plot rsq's measured on the in-fold data
    {
        if(!is.specified(col.infold.rsq))
            return()
        # recycle col.infold.rsq so can use different colors for different folds
        col.infold.rsq <- repl(col.infold.rsq, length(object$cv.list))
        for(ifold in seq_along(object$cv.list)) {
            infold.rsq <- object$cv.infold.rsq.tab[ifold,]
            if(jitter > 0)
                infold.rsq  <- jitter(infold.rsq, amount=jitter)
            scaled.rsq <- scale1(infold.rsq,  ylim[1], ylim[2])
            lines(scaled.rsq, col=col.infold.rsq[ifold], lty=1)
        }
    }
    draw.oof.rsqs <- function() # plot rsq's measured on the out-of-fold data
    {
        if(!is.specified(col.oof.rsq))
            return()
        # recycle col.oof.rsq so user can specify different colors for different folds
        col.oof.rsq <- repl(col.oof.rsq, length(object$cv.list))
        for(ifold in seq_along(object$cv.list)) {
            oof.rsq <- object$cv.oof.rsq.tab[ifold,]
            if(jitter > 0)
                oof.rsq  <- jitter(oof.rsq, amount=jitter)
            scaled.rsq <- scale1(oof.rsq,  ylim[1], ylim[2])
            lines(scaled.rsq, col=col.oof.rsq[ifold], lty=1)
        }
        if(is.specified(col.oof.labs)) {
            col.oof.labs <- repl(col.oof.labs, length(object$cv.list))
            x <- y <- labs <- NULL
            usr <- par("usr") # xmin, xmax, ymin, ymax
            for(ifold in seq_along(object$cv.list)) {
                oof.rsq <- object$cv.oof.rsq.tab[ifold,]
                oof.rsq <- oof.rsq[!is.na(oof.rsq)] # truncate NAs
                scaled.rsq <- scale1(oof.rsq,  ylim[1], ylim[2])
                y <- c(y, scaled.rsq[min(usr[2], length(oof.rsq))])
                x <- c(x, min(usr[2]-.1, length(oof.rsq)+.2))
                labs <- c(labs, substr(names(object$cv.list)[ifold], 5, 15))
            }
            cex <- .6
            text(x=x, y=TeachingDemos::spread.labs(y,
                            mindiff=1.2 * strheight("X")),
                 labels=labs, cex=cex, col=col.oof.labs[ifold], xpd=NA)
        }
        if(is.specified(col.pch.max.oof.rsq) || is.specified(col.pch.cv.rsq)) {
            for(ifold in seq_along(object$cv.list)) {
                oof.rsq <- object$cv.oof.rsq.tab[ifold,]
                scaled.rsq <- scale1(oof.rsq,  ylim[1], ylim[2])
                # show the max oof.rsq for this fold
                nterms <- which.max(oof.rsq)
                points(nterms, scale1(oof.rsq,  ylim[1], ylim[2])[nterms],
                       pch=1, col=col.pch.max.oof.rsq)
                # show the position of the cv.rsq's
                nterms <- length(object$cv.list[[ifold]]$selected.terms)
                points(nterms, scale1(oof.rsq,  ylim[1], ylim[2])[nterms],
                       pch=20, col=col.pch.cv.rsq, cex=.7)
            }
        }
    }
    draw.unused.preds <- function()  # plot nbr of used predictors
    {                               # nothing actually plotted if col.npreds=0
        nused.preds <- get.nused.preds.per.subset(object$dirs, object$prune.terms)
        nused.preds.vec <- scale1(nused.preds, 0, max.npreds)
        if(jitter > 0)  # 2*jitter seems to work better relative to jitter on GRSq
            nused.preds.vec <- jitter(nused.preds.vec, amount=2*jitter)
        else {
            # nudge max value to prevent overplot of maximum RSq(s)
            max <- max(nused.preds.vec)
            nused.preds.vec[nused.preds.vec == max] <- max + max / 150
        }
        lines(nused.preds.vec, type="l", col=col.npreds, lty=lty.npreds)
    }
    draw.vline.at.max.mean.oof.rsq <- function()
    {
        if(!is.specified(col.mean.oof.rsq) || !is.specified(col.oof.vline))
            return()
        x <- xnudge <- which.max(mean.oof.rsq.per.subset)
        # possibly nudge right to prevent overplot of grsq.line
        if(x == which.min(object$gcv.per.subset))
            xnudge <- xnudge + nterms.on.horiz.axis / 150
        # possibly nudge to prevent overplot of grid
        if(is.specified(grid.col))
            xnudge <- xnudge + nterms.on.horiz.axis / 150
        abline(v=xnudge, col=col.oof.vline, lty="12", lwd=1.2)
    }
    show.max.mean.oof.rsq <- function()
    {
        if(!is.specified(col.mean.oof.rsq) || !is.specified(col.oof.vline))
            return()
        x <- which.max(mean.oof.rsq.per.subset)
        if(which.min(object$gcv.per.subset) == x)
            return() # don't overplot (see show.max.grsq)
        usr <- par("usr")
        text.on.white(x, usr[3] + strheight("X"), x,
                      cex=.8, col=col.oof.vline, xmar=.05)
    }
    draw.mean.infold.rsq <- function()
    {
        if(!is.specified(col.mean.infold.rsq))
            return()
        lines(scale1(mean.infold.rsq.per.subset, ylim[1], ylim[2]),
              col=col.mean.infold.rsq, lwd=lwd)
    }
    draw.mean.oof.rsq <- function()
    {
        if(!is.specified(col.mean.oof.rsq))
            return()
        lines(scale1(mean.oof.rsq.per.subset, ylim[1], ylim[2]),
              col=col.mean.oof.rsq, lwd=lwd)
    }
    draw.rsq <- function()
    {
        if(jitter > 0)
            rsq.vec  <- jitter(rsq.vec, amount=jitter)
        lines(scale1(rsq.vec,  ylim[1], ylim[2]), col=col.rsq, lty=lty.rsq)
    }
    draw.grsq <- function()
    {
        if(jitter > 0)
            grsq.vec <- jitter(grsq.vec, amount=jitter)
        y <- scale1(grsq.vec, ylim[1], ylim[2])
        lines(y, col=col.grsq, lwd=lwd)
        # if pmethod=="cv", draw a circle at the selected model
        if(object$pmethod=="cv" && !is.null(mean.oof.rsq.per.subset)) {
            x <- length(object$selected.terms)
            points(x, y[x], col=col.grsq, lwd=lwd, pch=1)
        }
    }
    draw.vline.at.max.grsq <- function()
    {
        if(!is.specified(col.vline))
            return()
        x <- xnudge <- which.min(object$gcv.per.subset)
        # possibly nudge to prevent overplot of grid
        if(is.specified(grid.col))
            xnudge <- xnudge + nterms.on.horiz.axis / 150
        abline(v=xnudge, col=col.vline, lty=lty.vline, lwd=1.2)
        # possibly plot a colored marker at the top of the above line
        # (this is used by plot.earth.models when plotting multiple models)
        if(is.specified(col.vseg))
            points(x=xnudge, y=1.02, col=col.vseg, pch=6)
    }
    show.max.grsq <- function()
    {
        if(!is.specified(col.vline) || is.specified(col.vseg))
            return()
        x <- which.min(object$gcv.per.subset)
        usr <- par("usr")
        text.on.white(x, usr[3] + strheight("X"), x,
                      cex=.8, col=col.vline, xmar=.05)
    }
    draw.vline.at.max.nterms <- function()
    {
        if(object$pmethod != "none" || !is.specified(col.vline))
            return()
        # nrow(object$prune.terms) is nk
        x <- xnudge <- nrow(object$prune.terms)
        # possibly nudge to prevent overplot of grid
        if(is.specified(grid.col))
            xnudge <- xnudge + nterms.on.horiz.axis / 150
        # possibly nudge to prevent overplot of max.grsq line
        if(which.min(object$gcv.per.subset) == x)
            xnudge <- xnudge + nterms.on.horiz.axis / 150
        abline(v=xnudge, col=col.vline, lty=2, lwd=1.2)
    }
    show.max.nterms <- function()
    {
        if(object$pmethod != "none" || !is.specified(col.vline))
            return()
        x <- nrow(object$prune.terms)
        if(which.min(object$gcv.per.subset) == x ||
           (!is.null(mean.oof.rsq.per.subset) &&
                which.max(mean.oof.rsq.per.subset) == x)) {
            # don't overplot (see show.max.grsq and show.max.mean.oof.rsq)
            return()
        }
        usr <- par("usr")
        text.on.white(x, usr[3] + strheight("X"),
                      x, cex=.8, col=col.vline, , xmar=.05)
    }
    draw.legend <- function(...)
    {
        # return TRUE if "over" lines obscure "under" lines
        is.obscured <- function(under, over)
        {
            len <- min(length(under), length(over))
            under <- under[1:len]
            over  <- over[1:len]
            i <- under >= ylim[1] & under <= ylim[2]
            i[is.na(under) | is.na(over)] <- FALSE # ignore NAs
            nobscured <- sum(abs(under[i] - over[i]) < (ylim[2] - ylim[1]) / 100)
            nobscured > .8 * sum(i)
        }
        update.legend <- function(text, col=1, lty=1, lwd=1, vert=FALSE, pch=NA)
        {
            if(is.null(legend.text)) { # first time?
                if(text == "")         # spacer between entries?
                    return()           # ignore space when first entry
                legend.text <<- text   # note <<- not <-
                legend.col  <<- col
                legend.lty  <<- lty.as.char(lty)
                legend.lwd  <<- lwd
                legend.vert <<- vert
                legend.pch  <<- pch
            } else {
                legend.text <<- c(legend.text, text)
                legend.col  <<- c(legend.col, col)
                legend.lty  <<- c(legend.lty, lty.as.char(lty))
                legend.lwd  <<- c(legend.lwd, lwd)
                legend.vert <<- c(legend.vert, vert)
                legend.pch  <<- c(legend.pch, pch)
            }
        }
        #--- draw.legend starts here
        # The is.obscured code assumes that plot order is rsq, mean.oof.rsq, grsq.
        # Obscuring of or by infold.rsq is not yet handled.
        if(!is.null(legend.pos) && !is.specified(legend.pos))
            return()
        legend.text <- legend.col <- legend.lty <- legend.lwd <- NULL
        legend.vert <- legend.pch <- NULL
        full.model <- if(show.cv.data) " (full model)" else ""
        if(is.specified(col.vline) && object$pmethod == "none") {
            update.legend("selected model (pmethod=none)",
                col.vline, lty=2, lwd=1.2, vert=TRUE)
            update.legend("", 0) # dummy entry to leave a vertical space
        }
        if(is.specified(col.grsq))
            update.legend(paste0("GRSq", full.model), lwd=lwd)
        if(is.specified(col.vline))
            update.legend(
                if(object$pmethod != "none" && object$pmethod != "cv")
                    "selected model"
                else
                    "max GRSq",
                col.vline, lty.vline, lwd=1.2, vert=TRUE)
        if(is.specified(col.rsq)) {
            RSq.string <- if(show.cv.data) "RSq (full model)" else "RSq"
            if(is.specified(col.grsq) && is.obscured(rsq.vec, grsq.vec))
                text <- paste0(RSq.string, " (obscured)")
            else if(is.specified(col.mean.oof.rsq) &&
                    is.obscured(rsq.vec, mean.oof.rsq.per.subset))
                text <- paste0(RSq.string, " (obscured)")
            else
                text <- RSq.string
            update.legend(text, col.rsq, lty.rsq)
        }
        added.space <- FALSE
        # We draw the infold legend above the oof legend because the infold
        # curves are usually above the oof curves.
        if(is.specified(col.mean.infold.rsq)) {
            text <- "mean in-fold RSq"
            update.legend("", 0) # dummy entry to leave a vertical space
            added.space <- TRUE
            update.legend(text, col.mean.infold.rsq, lwd=lwd)
        }
        if(is.specified(col.infold.rsq)) {
            if(!added.space)
                update.legend("", 0) # dummy entry to leave a vertical space
            update.legend("in-fold RSq", col.infold.rsq[1])
        }
        if(is.specified(col.mean.oof.rsq)) {
            if(is.specified(col.grsq) && is.obscured(mean.oof.rsq.per.subset, grsq.vec))
                text <- "mean out-of-fold RSq (obscured)"
            else
                text <- "mean out-of-fold RSq"
            update.legend("", 0) # dummy entry to leave a vertical space
            added.space <- TRUE
            update.legend(text, col.mean.oof.rsq, lwd=lwd)
            if(is.specified(col.oof.vline))
                update.legend("max mean out-of-fold RSq", col.oof.vline, lty="12",
                              lwd=1.2, vert=TRUE)
        }
        if(is.specified(col.oof.rsq)) {
            if(!added.space)
                update.legend("", 0) # dummy entry to leave a vertical space
            update.legend("out-of-fold RSq", col.oof.rsq[1])
        }
        if(is.specified(col.npreds)) {
            if(added.space)
                update.legend("", 0) # dummy entry to leave a vertical space
            update.legend(paste0("nbr preds", full.model), col.npreds, lty.npreds)
        }
        if(is.specified(col.oof.vline) && object$pmethod == "cv")
            update.legend("selected model", col=col.grsq, lty=NA, lwd=lwd, pch=1)
        legend.cex <- get.earth.legend.cex(legend.text, ...)
        legend.inset <- 0
        if(is.null(legend.pos)) { # auto?
            if(max.nterms == 2) {
                legend.x <- "topleft"
                legend.inset <- c(.02, .02)
            } else {
                legend.x <- "bottomright"
                # legend y offset allows vertical lines and text to be visible below legend
                legend.inset <- c(.02, max(.05, 2 * strheight("X")))
            }
            legend.y <- NULL
            usr <- par("usr") # xmin, xmax, ymin, ymax
        } else { # user specified legend position
            legend.x <- legend.pos[1]
            legend.y <- NULL
            if(length(legend.pos) == 1) # presumably something like "topright"
                legend.inset <- c(.02, .02)
            else
                legend.y <- scale1(legend.pos[2], ylim[1], ylim[2])
        }
        if(max.nterms == 1)
            text.on.white(usr[1] + 2 * strwidth("X"),
                          usr[4] - 2 * strheight("X"),
                          "intercept-only model",
                          adj=0)
        else if(!is.null(object$nprune))
            text.on.white(usr[1] + .5 * strwidth("X"),
                          usr[4] - .6 * strheight("X"),
                          paste0("nprune ", object$nprune),
                          adj=0, cex=.8 * par("cex"))
        elegend(x=legend.x, y=legend.y, bg="white", legend=legend.text,
                col=legend.col, lty=legend.lty, lwd=legend.lwd,
                vert=legend.vert, pch=legend.pch,
                cex=legend.cex, xpd=NA, inset=legend.inset)
    }
    #--- earth_plotmodsel starts here ---
    object <- x
    remove(x) # prevent confusion with the x matrix
    main <- dot("main", ...)
    if(!is.specified(main))
        main <- if(NCOL(object$residuals) > 1) "Model Selection (all responses)"
                else                           "Model Selection"
    stopifnot.string(main, allow.empty=TRUE)
    if(is.null(object$prune.terms)) { # no prune data?
        if(!add)
            plot(c(0,1), col=0, xlab="", ylab="", main=main)
        legend(x=1, y=1, bty="n", legend=c("No model selection data", "",
               "Run update.earth() to generate", "model selection data"))
        return(NULL)
    }
    check.numeric.scalar(jitter)
    stopifnot(jitter >= 0)
    if(jitter > .1)
        stop0("'jitter' ", jitter , " is too big, try something like jitter=0.01")
    if(!is.specified(lty.grsq))   col.grsq <- 0
    if(!is.specified(lty.rsq))    col.rsq <- 0
    if(!is.specified(lty.npreds)) col.npreds <- 0
    if(!is.specified(lty.vline))  col.vline <- 0
    grid.col <- dot("grid.col col.sel.grid", ...)
    ylim <- get.model.selection.ylim(object,
                ylim=dot("ylim", DEF=NULL, ...), col.grsq=1, col.rsq,
                col.mean.oof.rsq, col.oof.rsq, col.mean.infold.rsq, col.infold.rsq)
    possibly.issue.cv.warning()
    if(is.null(object$cv.oof.rsq.tab)) # if no cv data available, force no display of cv data
        col.mean.oof.rsq <- col.oof.rsq <- col.mean.infold.rsq <- col.infold.rsq <- 0
    show.cv.data <-
        is.specified(col.mean.oof.rsq)  || is.specified(col.oof.rsq) ||
        is.specified(col.mean.infold.rsq) || is.specified(col.infold.rsq)
    show.non.cv.data <-
        is.specified(col.grsq) || is.specified(col.rsq) || is.specified(col.npreds)
    if(is.null(col.npreds)) # by default, show npreds if not show cv data
        col.npreds <- if(show.cv.data) 0 else 1
    rsq.vec  <- get.rsq(object$rss.per.subset, object$rss.per.subset[1])
    grsq.vec <- get.rsq(object$gcv.per.subset, object$gcv.per.subset[1])
    mean.oof.rsq.per.subset <- NULL
    if(is.specified(col.mean.oof.rsq))
        mean.oof.rsq.per.subset <-
            object$cv.oof.rsq.tab[nrow(object$cv.oof.rsq.tab),]
    if(is.specified(col.mean.infold.rsq))
        mean.infold.rsq.per.subset <-
            object$cv.infold.rsq.tab[nrow(object$cv.infold.rsq.tab),]
    lwd <- if(show.cv.data) 2 else 1 # want fat non-cv lines if plotting cv data
    nterms.on.horiz.axis <- max.nterms
    if(show.cv.data && !show.non.cv.data)
        nterms.on.horiz.axis <-
            min(nterms.on.horiz.axis, get.max.terms.of.fold.models(object))
    if(!add) {
        old.mar <- par("mar")
        if(is.specified(col.npreds) && old.mar[4] < 3.5) {
            # ensure right margin big enough for right axis
            on.exit(par(mar=old.mar))
            par(mar=c(old.mar[1:3], 3.5))
        }
        xlim <- get.model.selection.xlim(object, dot("xlim", ...),
                    mean.oof.rsq.per.subset, col.mean.oof.rsq, col.oof.vline)
        # set up so vertical scale is 0..1, horizontal is 0..nterms.on.horiz.axis
        plot(0:nterms.on.horiz.axis,
             (0:nterms.on.horiz.axis)/nterms.on.horiz.axis,
             type="n", main=main,
             xlab="Number of terms", ylab="", yaxt="n", xlim=xlim)
        left.axis()
        if(is.specified(col.npreds))
            right.axis()
        draw.selection.grid()
    }
    # note: if you change the plot order here, modify is.obscured code in draw.legend
    draw.infold.rsqs()
    draw.oof.rsqs()
    draw.unused.preds()
    draw.vline.at.max.grsq()
    draw.vline.at.max.mean.oof.rsq()
    draw.vline.at.max.nterms()
    draw.rsq()
    draw.mean.infold.rsq()
    draw.mean.oof.rsq()
    draw.grsq()
    show.max.grsq()
    show.max.mean.oof.rsq()
    show.max.nterms()
    draw.legend(...)
}
# Note: there is no string line type corresponding to 1, so this
# converts 1 to "1" which is an illegal lty, so must be specially
# handled in functions which use the lty string.

lty.as.char <- function(lty)
{
    char <- lty
    if(is.na(lty))
        char <- "NA"
    else if(is.numeric(lty)) {
        char <- NULL
        tab <- c("1", "44", "13", "1343", "73", "2262") # from par man page
        stopifnot(length(lty) > 0)
        for(i in seq_along(lty)) {
            stopifnot(lty[i] >= 1, lty[i] <= length(tab))
            char <- c(char, tab[lty[i]])
        }
    }
    char
}
get.earth.legend.cex <- function(legend.text, min.width=.4, min.cex=.4, ...)
{
    cex <- dot("legend.cex cex.legend", EX=c(0,1), NEW=1, ...)
    if(is.na(cex)) {
        longest.text <- legend.text[which.max(strwidth(legend.text))]
        longest.text <- paste0("AAAAAA ", longest.text) # incorporate line on left of legend
        # reduce cex until legend fits, but not more than min.cex
        cex <- .8
        while((width <- max(strwidth(longest.text, units="figure", cex=cex))) > min.width &&
                cex > min.cex)
            cex <- cex - .1
    }
    cex
}
get.model.selection.xlim <- function(object, xlim,
    mean.oof.rsq.per.subset, col.mean.oof.rsq, col.oof.vline)
{
    if(!is.specified(xlim)) { # not specified by the user?
        nk <- nrow(object$dirs)
        xmax <- 2 * which.min(object$gcv.per.subset) # nbr terms selected by GCV
        if(object$pmethod == "none")
            xmax <- max(xmax, nk)
        # if cross-validation vert line is plotted, include that too
        # following "if" matches that in draw.vline.at.max.mean.oof.rsq
        if(!is.null(mean.oof.rsq.per.subset) &&
                is.specified(col.mean.oof.rsq) && is.specified(col.oof.vline))
            xmax <- max(xmax, which.max(mean.oof.rsq.per.subset))
        xlim <- c(0, min(xmax + 3, nk))
    }
    xlim
}
# check ylim specified by user, and convert special values in ylim to actual vals

get.model.selection.ylim <- function(object, ylim, col.grsq, col.rsq,
                                     col.mean.oof.rsq=0, col.oof.rsq=0,
                                     col.mean.infold.rsq=0, col.infold.rsq=0)
{
    get.fold.min.max <- function()
    {
        min <- Inf
        max <- -Inf
        if(!is.null(object$cv.oof.rsq.tab) &&
                (is.specified(col.mean.oof.rsq) || is.specified(col.oof.rsq))) {
            # will be plotting oof.rsq, so must adjust axis limits for that
            min <- min(object$cv.oof.rsq.tab[,-1], na.rm=TRUE) # -1 to ignore intercept-only model
            max <- max(object$cv.oof.rsq.tab[,-1], na.rm=TRUE)
            # prevent outrageous axis scales caused by wayward cross-validation results
            max <- min(max, 2 * max(rsq))   # 2 is arb
            min <- max(min, -3)             # -3 is arb
        }
        if(!is.null(object$cv.infold.rsq.tab) &&
                (is.specified(col.mean.infold.rsq) || is.specified(col.infold.rsq))) {
            min <- min(min, object$cv.infold.rsq.tab[,-1], na.rm=TRUE)
            max <- max(max, object$cv.infold.rsq.tab[,-1], na.rm=TRUE)
            max <- min(max, 2 * max(rsq))
            min <- max(min, -3)
        }
        list(min=min, max=max)
    }
    #--- get.model.selection.ylim starts here ---
    if(is.null(ylim))
        ylim <- c(-1, -1)
    if(length(ylim) != 2)
        stop0("length(ylim) != 2")
    if(ylim[2] <= ylim[1] && ylim[2] != -1)
        stop0("ylim[2] <= ylim[1]")
    if(ylim[1] < -1 || ylim[1] >  1 || ylim[2] < -1 || ylim[2] >  1)
        stop0(paste0(
              "illegal ylim=c(", ylim[1], ",", ylim[2], ") in the earth model selection plot\n",
              "Allowed settings are from -1 to 1, with special values:\n",
              "  ylim[1] = -1 means use min(RSq, GRSq)\n",
              "  ylim[2] = -1 means use max(RSq, GRSq)\n"))
    if(ylim[1] == -1 || ylim[2] == -1) {
        grsq <- NULL
        if(is.specified(col.grsq))
            grsq <- get.rsq(object$gcv.per.subset, object$gcv.per.subset[1])
        rsq <- get.rsq(object$rss.per.subset, object$rss.per.subset[1])
        temp <- get.fold.min.max()
        if(!is.specified(col.rsq))
            rsq <- NULL
        if(ylim[1] == -1) {
            ylim[1] <- min(grsq[-1], rsq[-1], temp$min, na.rm=TRUE)
            # small model, treat specially so user sees context
            if(length(object$rss.per.subset) <= 3)
                ylim[1] <- min(0, ylim[1])
            ylim[1] <- max(-1, ylim[1]) # clamp minimum ylim at -1
        }
        if(ylim[2] == -1)
            ylim[2] <- max(grsq, rsq, temp$max, na.rm=TRUE)
    }
    # following code gives a decent y axis even with an intercept-only model
    if(abs(ylim[1] - ylim[2]) < 1e-6)
        ylim[2] <- ylim[1] + 1
    ylim
}
get.max.terms.of.fold.models <- function(object)
{
    tab <- object$cv.oof.rsq.tab
    stopifnot(!is.null(tab))
    stopifnot(nrow(tab) > 1)
    max.terms <- 0
    for(i in 1:(nrow(tab)-1)) # -1 to skip last (summary) row
        max.terms <- max(max.terms, sum(!is.na(tab[i,])))
    max.terms
}
# Given a list of objects, return a vector of strings.  Each string shows where
# the $call argument of the object differs from the call of the first object.
# (Thus the first object is used as a reference).

get.arg.strings <- function(
        objects,    # list of objects with $call arguments
        maxchars=16)
{
    # the gsub discards white space and the final ")"
    get.call <- function(iobj)
        gsub("[ \t\n)]",  "",  paste.collapse(format(objects[[iobj]]$call)))

    stopifnot(length(objects) > 0)
    call <- get.call(1)
    if(length(objects) == 1)
        return(substr(call, 1, maxchars))
    call2 <- get.call(2)
    i <- first.non.matching.arg(call, call2)
    if(i == 0)
        rval <- c("", "")
    else
        rval <- c(substr(call, i, i+maxchars), substr(call2, i, i+maxchars))
    if(length(objects) > 2)
        for(iobj in 3:(length(objects))) {
            call2 <- get.call(iobj)
            i <- first.non.matching.arg(call, call2)
            rval <- c(rval, if(i==0) "" else substr(call2, i, i+maxchars))
        }
    rval
}
# Return the position of the first non matching arg between two function call strings.
#
# More precisely, find the first non matching characters in s1 and s2.
# When it is found, step back until a comma or "(" is found.
# Return the index of the character after the comma or "(".
#
# Example: s1 = lm(formula=O3~.,data=ozone
#          s2 = lm(formula=O3~.-wind,data=ozone
#
#          return index of "formula=O3~.-wind,data=ozone"
#          because formula is the first argument with a differing argument

first.non.matching.arg <- function(s1, s2)
{
    len <- min(nchar(s1), nchar(s2))
    if(len == 0)
        return(0)
    for(i in 1:len)
        if(substr(s1, i, i) != substr(s2, i, i))
            break
    if(i == len || i == 1)  # no difference or all different?
        return(1)
    while(i >= 1 && substr(s2, i, i) != "," && substr(s2, i, i) != "(")
        i <- i - 1  # move backwards to character following comma or "("
    return(i+1)
}
