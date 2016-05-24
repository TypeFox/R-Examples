# plot.gbmx.R: plot gbm models
#
# This code is based on gbm 2.1.1 (April 2015).
#
# TODO maybe add arg to rescale errs e.g. RSquared rather than Squared Error
# TODO add right hand axis for OOB, or scale OOB to same units when possible
# TODO when selecting best n.trees, why is OOB smoothed but not test or CV?
# TODO maybe use dots package (and then won't need explicit xlim, ylim, etc. args)
# TODO add test of args like col and legend.x to plotmo test suite
# TODO if gbm calculated CV stddev across folds we could plot conf bands

plot.gbmx <- function(object=stop("no 'object' argument"),
    xlim = NULL,            # NULL means default
    ylim = NULL,
    col = c(1, 2, 3, 4),    # colors for train, test, CV, OOB
    legend.x = "topright",  # use legend.x=NA for no legend
    legend.y = NULL,
    legend.cex = .8,
    n.trees = NA,           # draw a vertical gray line at n.trees (for plotres)
    ...)                    # dots passed to plot.default
{
    ylab <- if (object$distribution$name =="pairwise")
                switch(object$distribution$metric,
                        conc="Fraction of concordant pairs",
                        ndcg="Normalized discounted cumulative gain",
                        map ="Mean average precision",
                        mrr ="Mean reciprocal rank",
                        stop0("bad object$distribution$metric"))
            else # not pairwise
                switch(substring(object$distribution$name,1,2),
                        ga="Squared error loss",
                        be="Bernoulli deviance",
                        po="Poisson deviance",
                        ad="AdaBoost exponential bound",
                        co="Cox partial deviance",
                        la="Absolute loss",
                        qu="Quantile loss",
                        mu="Multinomial deviance",
                        td="t-distribution deviance",
                        stop0("bad object$distribution$name"))

    if(!is.specified(xlim))
        xlim <- c(0, object$n.trees)
    if(!is.specified(ylim))
         ylim <- range(object$train.error, object$valid.error, object$cv.error,
                       na.rm=TRUE)
    col <- rep_len(col, 4) # recycle col if necessary

    plot(1:object$n.trees, object$train.error, type="n",
         xlab="Number of trees", ylab=ylab, xlim=xlim, ylim=ylim, ...)

    # draw n.trees vertical gray line first, so other plots go on top of it
    col.n.trees <- "darkgray"
    if(is.specified(n.trees))
        vert.line(n.trees, col.n.trees, 1, 0)

    legend.text <- legend.col <- legend.lty <- legend.vert <- legend.n <- NULL
    offset <- 0 # offset is to distinguish overplotted vertical lines

    # train
    lines(object$train.error, type="l", col=col[1])
    legend.text <- c(legend.text,
                     if(object$train.fraction == 1) "train"
                     else sprintf("train (frac %g)", object$train.fraction))
    legend.col  <- c(legend.col, col[1])
    legend.lty  <- c(legend.lty, 1)
    legend.vert <- c(legend.vert, FALSE)
    legend.n    <- which.min(object$train.error)

    # test (aka valid.error)
    if(object$train.fraction != 1) {
        min <- which.min(object$valid.error)
        vert.line(min, 2, 3, offset)
        offset <- offset + 1
        lines(object$valid.error, col=col[2])
        legend.text <- c(legend.text,
                         sprintf("test (frac %g)", 1-object$train.fraction))
        legend.col  <- c(legend.col, col[2])
        legend.lty  <- c(legend.lty, 1)
        legend.vert <- c(legend.vert, FALSE)
        legend.n    <- c(legend.n, min)
    }
    # CV
    if(!is.null(object$cv.error)) {
        stopifnot(object$cv.folds > 1)
        min <- which.min(object$cv.error)
        vert.line(min, col[3], 3, offset)
        offset <- offset + 1
        lines(object$cv.error, col=col[3])
        legend.text <- c(legend.text, sprintf("CV (%g fold)", object$cv.folds))
        legend.col  <- c(legend.col, col[3])
        legend.lty  <- c(legend.lty, 1)
        legend.vert <- c(legend.vert, FALSE)
        legend.n    <- c(legend.n, min)
    }
    # OOB
    if(object$bag.fraction != 1) {
        min <- draw.oob(object, offset, col[4])
        offset <- offset + 1
        legend.text <- c(legend.text,
                         if(is.na(min)) "OOB not plotted" else "OOB (rescaled)")
        legend.col  <- c(legend.col, col[4])
        legend.lty  <- c(legend.lty, 2)
        legend.vert <- c(legend.vert, FALSE)
        legend.n    <- c(legend.n, min)
    }
    if(is.specified(n.trees)) {
        legend.text <- c(legend.text, "predict n.trees")
        legend.col  <- c(legend.col, col.n.trees)
        legend.lty  <- c(legend.lty, 1)
        legend.vert <- c(legend.vert, TRUE)
        legend.n    <- c(legend.n, n.trees)
    }
    box() # replot box because vert.line overplots it slightly
    if(is.specified(legend.x))
        elegend(x=legend.x, y=legend.y,
                legend=legend.text, col=legend.col, lty=legend.lty,
                vert=legend.vert, # vert is supported by elegend but not by legend
                bg="white", cex=legend.cex)

    toplabs(legend.n, legend.text, legend.col, legend.cex)
}
vert.line <- function(x, col=1, lty=1, offset=0) # draw a vertical line at x
{
    usr <- par("usr") # xmin, xmax, ymin, ymax
    range <- usr[4] - usr[3]
    # small vertical offset (only visible for non-solid linetypes)
    offset <- 0.005 * offset * range
    lines(x=c(x, x), y=c(usr[3], usr[4]) - offset, col=col, lty=lty)
    lines(x=c(x, x), y=c(usr[3], usr[3]+.02*range), col=col, lty=1) # tick
}
draw.oob <- function(object, offset, col)
{
    if(all(!is.finite(object$oobag.improve))) { # TODO can this actually happen?
        warning0("plot.gbmx: cannot plot OOB curve (bad oobag.improve)")
        return(NA)
    }
    stopifnot(length(object$oobag.improve) == object$n.trees) # paranoia

    # to calculate min use same smoothing as gbm.perf, for compatibility
    # (but plot itself is unsmoothed, else tend to smooth away left part of curve)
    min <- min.smoothed.oob(object)
    vert.line(min, col, 3, offset)

    y <- -cumsum(object$oobag.improve)
    # scale y to fit in the existing plot (hack)
    # we scale so y takes up 80% of ylim (i.e. a 10% space at the top and bottom)
    y <- (y  - min(y)) / (max(y) - min(y))
    usr <- par("usr") # xmin, xmax, ymin, ymax
    yrange <- usr[4] - usr[3]
    y <- usr[3] + .1 * yrange + .8 * yrange * y

    lines(1:object$n.trees, y, col=col, lty=2)
    min
}
min.smoothed.oob <- function(object) # same algorithm as gbm.perf (gbm version 2.1.1)
{
  x <- 1:object$n.trees
  smoother <- loess(object$oobag.improve ~ x,
                    na.action=na.omit, # paranoia, prevent warnings from loess
                    enp.target=min(max(4, length(x) /10 ), 50 ))
  x[which.min(-cumsum(smoother$fitted))]
}
# print the best n for each plot along the top of the graph
toplabs <- function(legend.n, legend.text, legend.col, legend.cex)
{
    # don't print n for the training curve
    stopifnot(substring(legend.text[1], 1, 5) == "train")
    legend.col[1] <- 0

    # darker than darkgray "#A9A9A9" seems needed for text to be perceived as darkgray
    legend.col[legend.col == "darkgray"] <- "#707070"

    usr <- par("usr") # xmin, xmax, ymin, ymax

    # TODO spread.labs is buggy for horizontal labels (too much space)?
    text(x=TeachingDemos::spread.labs(legend.n,
                mindiff=legend.cex * max(strwidth(paste(legend.n))),
                min=0, max=usr[2]),
         y=usr[4] + .03 * (usr[4] - usr[3]), # just above top of plot
         labels=legend.n, col=legend.col, cex=legend.cex, xpd=NA)
}
