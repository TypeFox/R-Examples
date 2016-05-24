# plotd.R: plot densities of class conditional predicted values
#
# TODO: allow newdata so can plot not only with the training data
# TODO: allow freq arg for histograms

plotd <- function(object,   # object is a model object
    hist   = FALSE,         # FALSE to use density(), TRUE to use hist()
    type   = NULL,          # NULL gets changed to a value which is passed on to predict
    nresponse = NULL,       # which response, for multiple response models, NULL for all
    dichot = FALSE,
    trace  = FALSE,
    xlim   = NULL,          # NULL means auto
    ylim   = NULL,          # NULL means auto
    jitter = FALSE,
    main   = NULL,          # graph caption
                            #   "string"  string
                            #   ""        no caption
                            #   NULL      generate a caption from x$call
    xlab = "Predicted Value",
    ylab = if(hist) "Count" else "Density",
    lty  = 1,               # linetypes for the plotted lines
    col  = c("gray70", 1, "lightblue", "brown", "pink", 2, 3, 4), # cols for plotted lines
    fill = if(hist) col[1] else 0, # fill color for first hist/density plot
    breaks = "Sturges",     # following passed on to hist, only used if hist=TRUE
    labels = FALSE,
    kernel = "gaussian",   # following passed on to density, only used if hist=FALSE
    adjust = 1,
    zero.line = FALSE,
    legend = TRUE,          # TRUE to draw a legend
    legend.names = NULL,    # NULL means auto, else specify a vector of strings,
    legend.pos = NULL,      # NULL means auto, else specify c(x,y) in user coords
    cex.legend = .8,        # cex for legend
    legend.bg = "white",    # bg for legend
    legend.extra = FALSE,   # print number in each class in legend
    vline.col = 0,          # color of vertical line, use NULL or 0 for no line
    vline.thresh = .5,      # horizontal position of vertical line
    vline.lty = 1,          # lty of vertical line
    vline.lwd = 1,          # lwd of vertical line
    err.thresh=vline.thresh, # thresh for "error areas"
    err.col=0,              # col shading of "error areas"
    err.border=0,
    err.lwd=1,
    xaxt = "s",
    yaxt = "s",
    xaxis.cex = 1,
    sd.thresh = 0.01,
    ...)                    # passed to predict
{
    init.global.data()
    on.exit(init.global.data()) # release memory on exit
    object.name <- short.deparse(substitute(object))
    trace <- as.numeric(check.integer.scalar(trace, logical.ok=TRUE))
    # Associate the model environment with the object.
    # (This is instead of passing it as an argument to plotmo's data access
    # functions.  It saves a few hundred references to model.env in the code.)
    attr(object, ".Environment") <- get.model.env(object, object.name, trace)
    temp <- plotmo::plotmo_prolog(object, object.name, trace, ...)
        object  <- temp$object
        my.call <- temp$my.call
    hist         <- check.boolean(hist)
    # dicho      <- check.boolean(dichot) # can't use because we use missing(dichot) below
    trace        <- as.numeric(check.integer.scalar(trace, logical.ok=TRUE))
    jitter       <- check.boolean(jitter)
    labels       <- check.boolean(labels)
    zero.line    <- check.boolean(zero.line)
    legend       <- check.boolean(legend)
    legend.extra <- check.boolean(legend.extra)

    type <- plotmo::plotmo_type(object, trace, "plotmo", type, ...)

    yhat.per.class <- get.yhat.per.class(object, object.name,
                            type, nresponse, dichot, trace, ...)

    nclasses <- length(yhat.per.class)

    # get densities

    densities <- NULL
    for(iclass in seq_len(nclasses))
        if(!hist)
            densities[[iclass]] <- density(yhat.per.class[[iclass]],
                                           kernel=kernel, adjust=adjust)
        else {
            densities[[iclass]] <- hist(yhat.per.class[[iclass]],
                                        breaks=breaks, plot=FALSE)
            # need x and y components so hist can be treated uniformly with density
            densities[[iclass]]$x <- densities[[iclass]]$breaks
            densities[[iclass]]$y <- densities[[iclass]]$counts
        }

    # get x limits of plot

    if(is.null(xlim)) {
        min1 <- Inf
        max1 <- -Inf
        for(iclass in seq_len(nclasses)) {
             min1 <- min(min1, densities[[iclass]]$x)
             max1 <- max(max1, densities[[iclass]]$x)
        }
        xlim <- c(min1, max1)
    }
    if(length(xlim) != 2)
        stop0("length(xlim) != 2")
    xspan <- xlim[2] - xlim[1]

    # sanity check the ranges of each class, issue warnings if need be

    degenerate <- logical(nclasses)
    for(iclass in seq_len(nclasses)) {
        range <- range(yhat.per.class[[iclass]])
        if(sd(yhat.per.class[[iclass]]) < sd.thresh) {
            warning0("standard deviation of '", names(yhat.per.class)[iclass],
                     "' density is ", sd(yhat.per.class[[iclass]]),
                     ",  density is degenerate?")
            degenerate[iclass] <- TRUE
        }
    }
    # add jitter and get ymax for plot

    ymax <- 1e-6
    if(is.logical(jitter) && jitter)
        jitter <- xspan / 100
    for(iclass in seq_len(nclasses)) {
        if(jitter) {
            if(hist)
                densities[[iclass]]$breaks <-
                    densities[[iclass]]$breaks + iclass * jitter
            else
                densities[[iclass]]$x <- densities[[iclass]]$x + iclass * jitter
        }
        ymax <- max(ymax, densities[[iclass]]$y)
    }
    if((is.logical(labels) && labels) || is.character(labels))
        ymax <- 1.1 * ymax # hack to make space for labels
    if(!is.null(ylim)) {
        if(length(ylim) != 2)
            stop0("length(ylim) != 2")
        ymax <- ylim[2]
        if(ymax <= 0)
            stop0("ylim[2] <= 0")
        if(ylim[1] != 0)
            warning0("ignoring ylim[1], treating it as 0")
    }
    # expand lty and other arguments if necessary

    if(length(lty) < nclasses)
        lty <- repl(lty, nclasses)
    if(length(col) < nclasses)
        col <- repl(col, nclasses)

    if(is.null(main)) { # auto generate main?
        main <- paste0(object.name, " ", paste.collapse(type),
                      if(missing(nresponse)) ""
                      else paste0(" nresp=", paste(nresponse, collapse=",")),
                      if(missing(dichot)) ""
                      else paste0(" dichot=", dichot))
        main <- paste.trunc(main)
    }
    if(!is.null(my.call)) {
        main <- paste0(main, "\n", paste.trunc(my.call))
        old.cex.main <- par("cex.main")
        on.exit(par(cex.main=old.cex.main), add=TRUE)
        par(cex.main=1)
    }
    # we draw our own x axis for type="class"
    # xlims are wrong for histograms if a density is degenerate hence the test
    # TODO weird behaviour in hist? hist gives a 0 lower x val if density degenerate
    draw.own.axis <- hist && xaxt != "n" && !is.na(pmatch(type, "class")) &&
                     !any(degenerate)
    if(draw.own.axis)
        xaxt <- "n"

    # plot the first graph

    ifirst <- 1 # index of first non-degenerate class, 1 if all degenerate
    for(iclass in seq_len(nclasses))
        if(!degenerate[iclass]) {
            ifirst <- iclass
            break
        }
    if(hist) { # plot.histogram
        plot(densities[[ifirst]], xlim=xlim, ylim=c(0, ymax),
             main=main, xlab=xlab, ylab=ylab, lty=lty[ifirst],
             border=col[ifirst], xaxt=xaxt, yaxt=yaxt,
             col=if(ifirst==1) fill else 0) # fill color
        draw.labels(densities[[ifirst]], labels, cex.legend)
    } else {   # plot.density
        plot(densities[[ifirst]], xlim=xlim, ylim=c(0, ymax), col=col[ifirst],
             main=main, xlab=xlab, ylab=ylab, lty=lty[ifirst],
             zero.line=zero.line, xaxt=xaxt, yaxt=yaxt)
        if(ifirst == 1 && is.specified(fill))
            polygon(densities[[1]], col=fill, border=col[1])
    }
    if(draw.own.axis) {
        at <- xlim[1]:xlim[2]
        # hack to adjust leftmost label TODO not quite right, but ok
        at[1] <- at[1] + strwidth(names(yhat.per.class)[1], "user")
        mtext(names(yhat.per.class), side=1, at=at, font=2, adj=1, cex=xaxis.cex)
    }
    # optional error region shading

    if(!hist && any1(err.col))
        draw.err.col(densities, err.thresh, err.col, err.border, err.lwd)

    # overlay the graphs

    for(iclass in seq_len(nclasses))
        if(!degenerate[iclass]) {
            if(!hist)  # lines.density
                lines(densities[[iclass]], col=col[iclass], lty=lty[iclass])
            else {     # lines.histogram
                lines(densities[[iclass]], col=NULL,
                      border=col[iclass], lty=lty[iclass])
                draw.labels(densities[[iclass]], labels, cex.legend)
            }
        }
    # optional vertical line at vline.thresh

    if(!is.null(vline.col) && is.specified(vline.col))
        abline(v=vline.thresh, col=vline.col, lty=vline.lty, lwd=vline.lwd)

    # Redo optional error region shading if it has borders, because the
    # borders go on top of the other plotted lines.

    if(any1(err.border))
        draw.err.col(densities, err.thresh, err.col, err.border, err.lwd)

    # optional legend
    if(legend)
        draw.legend(densities, degenerate, yhat.per.class, ymax,
                    hist, xlim, col, fill, lty, legend.names,
                    legend.pos, cex.legend, legend.bg, legend.extra)
    invisible(yhat.per.class)
}
# add histogram labels --- lifted from plot.histogram and
# tweaked to (i) use cex and (ii) not draw zero counts
draw.labels <- function(x, labels, cex)
{
    if((is.logical <- is.logical(labels) && labels) || is.character(labels)) {
        stopifnot(!is.null(x$counts)) # plotd hist supports only counts, not densities
        if(is.logical) {
            labels <- format(x$counts)
            labels[x$counts == 0] <- ""
        }
        text(x$mids, x$counts, labels=labels, adj=c(.5, -.5), cex=cex)
    }
}
is.numlog <- function(x)
    is.numeric(x) || is.logical(x)

is.lda.or.qda <- function(object) # allows hacks for lda and qda specific code
    inherits(object, "lda") || inherits(object, "qda")

# return the predictions for each class, in a list with each class named
# nomeclature: y is the observed response, yhat is the predicted response

get.yhat.per.class <- function(object, object.name, type, nresponse, dichot, trace, ...)
{
    temp <- get.plotd.data(object, type, nresponse, trace, ...)
        y             <- temp$y
        yhat          <- temp$yhat
        colnames.yhat <- temp$colnames.yhat
        nresponse     <- temp$nresponse

    yhat1 <- if(length(dim(yhat)) > 1) yhat[,1] else yhat # TODO could probably delete
    yhat.per.class <- list() # will put per-class predicted vals in here
    ylevs <- levels(y) # null if y is not a factor
    nlevs <- 0
    if(NCOL(yhat) == 1) {
        #---single column yhat--------------------------------------------------
        trace2(trace, "single column yhat\n")
        if(is.factor(y) && (is.numlog(yhat1) || is.factor(yhat1))) {
            nlevs <- nlevels(y)
            if(!is.numlog(yhat1) && !is.factor(yhat1))
                cannot.plot.this.response(y, yhat, colnames.yhat, type, nlevs)
            stopifnot(length(ylevs) > 1)
            if(!is.na(pmatch(type, "class")))
                dichot <- FALSE # no dichot for type="class"
            if(length(ylevs) == 2 || dichot) {
                ylev1 <- ylevs[1]
                if(length(ylevs) == 2) {
                    other.level.name <- paste(ylevs[2])
                    observed.string <- "two-level factor"
                } else {
                    other.level.name <- paste("not", ylev1)
                    observed.string <- "multi-level factor, but dichot"
                }
                trace.response.type(trace, type,
                    observed=observed.string,
                    predicted="numeric or logical vector",
                    "CLASS1 predicted[observed == ", ylev1,
                    "], CLASS2 predicted[observed == ", other.level.name, "]")
                yhat.per.class[[1]] <- yhat1[y == ylev1]
                check.min(yhat.per.class[[1]], ylev1)
                yhat.per.class[[2]] <- yhat1[y != ylev1]
                check.min(yhat.per.class[[2]], other.level.name)
                names(yhat.per.class) <- get.prefixed.names(c(ylev1, other.level.name),
                                                            yhat, colnames.yhat, nresponse)
            } else {
                trace.response.type(trace, type,
                    observed="multi-level factor",
                    predicted="numeric or logical vector",
                    "predicted[observed == level] for ",
                    length(ylevs), " levels")
                for(iclass in seq_along(ylevs)) {
                    lev <- ylevs[iclass]
                    yhat.per.class[[iclass]] <- yhat1[y == lev]
                    check.min(yhat.per.class[[iclass]], lev)
                }
                names(yhat.per.class) <- get.prefixed.names(ylevs, yhat, colnames.yhat, nresponse)
            }
        } else if(NCOL(y) == 2 && is.numlog(y[,1]) && is.numlog(y[,2])) {
            trace.response.type(trace, type,
                observed="numeric or logical vector",
                predicted="two-column numeric",
                "CLASS1 observed[,1] <= observed[,2], ",
                "CLASS2 observed[,1] > observed[,2]")
            # split into two classes based on relative sizes of columns of y
            yhat.per.class[[1]] <- yhat1[y[,1] <= y[,2]]
            check.min(yhat.per.class[[1]], "observed[,1] <= observed[,2]")
            yhat.per.class[[2]] <- yhat1[y[,1] > y[,2]]
            check.min(yhat.per.class[[2]], "observed[,1] > observed[,2]")
            names(yhat.per.class) <-
                get.binary.class.names(yhat, colnames.yhat, object$fitted.values, c("FALSE", "TRUE"))
        } else if(NCOL(y) == 1 && is.numlog(y)) {
            th <- get.thresh(y, "response")
            trace.response.type(trace, type,
                observed="numeric or logical vector",
                predicted="numeric or logical vector",
                "CLASS1 ", th$text.le, ", CLASS2 ", th$text.gt)
            yhat.per.class[[1]] <- yhat1[y <= th$thresh]
            check.min(yhat.per.class[[1]], th$text.le)
            yhat.per.class[[2]] <- yhat1[y > th$thresh]
            check.min(yhat.per.class[[2]], th$text.gt)
            names(yhat.per.class) <- NULL
            if(th$thresh == 0)
                names(yhat.per.class) <-
                    get.binary.class.names(yhat, colnames.yhat, object$fitted.values,
                                           c(th$text.le, th$text.gt))
            else
                names(yhat.per.class) <- c(th$text.le, th$text.gt)
        } else
            cannot.plot.this.response(y, yhat, colnames.yhat, type, nlevs)
  } else {
        #---multiple column yhat------------------------------------------------
        trace2(trace, "multiple column yhat\n")
        if(!is.numeric(yhat[,1]))
            cannot.plot.this.response(y, yhat, colnames.yhat, type, nlevs)
        if(NCOL(y) == 1 && is.null(ylevs))
            ylevs <- as.numeric(names(table(y))) # use numeric levels like a factor
        nlevs <- length(ylevs)
        if(NCOL(y) == 1 && nlevs == NCOL(yhat)) {
            if(is.factor(y))
                trace.response.type(trace, type,
                    observed="factor",
                    predicted=
"multicolumn numeric, ncol(predicted) == nlevels(observed)",
                    "observed==level for each level in observed response")
            else
                trace.response.type(trace, type,
                    observed="factor",
                    predicted=
"multicolumn numeric, ncol(predicted) == nbr.of.unique.vals.in.observed",
                    "observed==val for each unique val in observed response")
            for(iclass in seq_len(ncol(yhat))) {
                lev <- ylevs[iclass]
                yhat.per.class[[iclass]] <- yhat[y == lev, iclass]
                check.min(yhat.per.class[[iclass]], lev)
                if(length(yhat.per.class[[iclass]]) == length(yhat[,iclass]))
                    stop0("no occurrences of ", lev,
                          " in the observed response")
            }
            if(!is.null(colnames.yhat))
                names(yhat.per.class) <- colnames.yhat
            else
                names(yhat.per.class) <- ylevs
#         } else if(NCOL(y) == 1) { # nlevs != NCOL(yhat))
#           trace.response.type(trace, type,
#                   observed="factor",
#                   predicted="multicolumn numeric; ncol(predicted) != nlevels(observed)",
#                   "each column of predicted response is a group")
#             for(iclass in seq_len(ncol(yhat)))
#                 yhat.per.class[[iclass]] <- yhat[, iclass]
#             if(!is.null(colnames.yhat))
#                 names(yhat.per.class) <- colnames.yhat
#             else
#                 names(yhat.per.class) <-  paste0(type, "[,", seq_len(ncol(yhat)), "]")
        } else if(is.numeric(y) && NCOL(y) == NCOL(yhat)) {
            th <- get.thresh(y, "response")
            trace.response.type(trace, type,
                observed="multicolumn numeric",
                predicted=
"multicolumn numeric with same number of columns as observed response",
                th$text.gt, "for each column of observed response")
            for(iclass in seq_len(ncol(yhat))) {
                yhat.per.class[[iclass]] <- yhat[y[,iclass] > th$thresh, iclass]
                check.min(yhat.per.class[[iclass]], th$text.gt)
                if(length(yhat.per.class[[iclass]]) == length(yhat[,iclass]))
                    stop0("no occurrences of ", th$text.le,
                          " in the observed response")
            }
            names(yhat.per.class) <- get.class.names(y, yhat, colnames.yhat, object$fitted.values)
        } else
            cannot.plot.this.response(y, yhat, colnames.yhat, type, nlevs,
                "Remedy: use the \"nresponse\" argument to select ",
                "just one column of the predicted response\n")
    }
    nchar <- max(nchar(names(yhat.per.class)))
    for(iclass in seq_along(yhat.per.class)) {
        # density needs numeric
        yhat.per.class[[iclass]] <- as.numeric(yhat.per.class[[iclass]])
        if(trace >= 1) {
            trace2(trace, "\n")
            print_summary(yhat.per.class[[iclass]],
                          sprintf("predicted.response.per.class[%-*s]",
                                  nchar, names(yhat.per.class)[iclass]),
                          trace=max(2, trace), details=if(trace>=2) 2 else 0)
        }
    }
    if(trace >= 1)
        cat("\n")
    yhat.per.class
}
get.plotd.data <- function(object, type, nresponse, trace, ...)
{
    # TODO this routine is bit messy
    #      if it were unified with plotmo_meta we would support more models

    # assignInMyNamespace("trace.call.global", trace)
    y <- get.observed.response(object)
    if(trace >= 2) {
        print_summary(y, "observed response", trace=2)
        trace2(trace, "\n")
    }
    yhat <- plotmo::plotmo_predict(object, newdata=NULL, nresponse=NULL, type,
                                   expected.levs=NULL, trace, ...)$yhat
    # assignInMyNamespace("trace.call.global", 0)
    if(is.character(yhat[,1])) {
        if(trace >= 1)
            printf("convert character yhat to factor\n")
        expected.levs <- plotmo::plotmo_resplevs(object, NULL, y, trace)
        yhat[,1] <- factor(yhat[,1], levels=expected.levs)
    }
    colnames.yhat <- colnames(yhat)
    if(!is.null(nresponse)) {
        nresponse <- plotmo::plotmo_nresponse(yhat, object, nresponse, trace,
                                        sprintf("predict.%s", class(object)[1]),
                                        type)
        if(NCOL(yhat) > 1) {
            yhat <- yhat[, nresponse]
            if(is.data.frame(yhat)) # TODO needed for fda type="hier", why?
                yhat <- as.matrix(yhat)
             print_summary(yhat,
                paste("predict after selecting nresponse", nresponse), trace)
            trace2(trace, "\n")
        }
    }
    list(y             = y,
         yhat          = yhat,
         colnames.yhat = colnames.yhat,
         nresponse     = nresponse)
}
# get the original observed response (it's needed to determine correct classes)

get.observed.response <- function(object)
{
    if(!is.null(object$call$formula)) {
        # get y from formula and data used in original call
        data <- get.update.arg(NULL, "data", object, parent.frame(), FALSE)
        call <- object$call
        m <- match(c("formula", "data", "na.action", "subset"), names(call), 0)
        mf <- call[c(1, m)]
        mf[[1]] <- as.name("model.frame")
        mf <- eval(mf, model.env(object))
        y <- model.response(mf, "any")  # "any" means factors are allowed
        if(NCOL(y) == 1 && is.numlog(y)) {
            # turn into a matrix so we have the column name
            names(y) <- NULL # we don't need row names
            y <- as.matrix(y)
            colnames(y) <- colnames(mf)[attr(object$terms,"response")]
        }
    } else if(is.lda.or.qda(object)) { # hack for lda and qda, get grouping arg
        y <- eval(object$call[[3]], model.env(object))
        # sanity check
        if(NCOL(y) != 1 || length(y) < 3 || (!is.numeric(y) && !is.factor(y)))
            stop0("cannot get \"grouping\" argument from object$call")
    } else
        y <- get.update.arg(NULL, "y", object, parent.frame(),
                            trace1=FALSE, reeval=FALSE)

    if(is.lda.or.qda(object))
        y <- as.factor(y)  # to make plotd handle response appropriately
    y
}
check.min <- function(x, ...)
{
    len <- length(x)
    if(len == 0)
        warning0("no occurrences of ", paste0(...),
                 " in the observed response")
    else if(len < 3) # 3 is arbitrary
        warning0("only ", len, " occurrences of ", paste0(...),
                 " in the observed response")
}
get.binary.class.names <- function(yhat, colnames.yhat, fitted.values, last.resort)
{
    if(!is.null(colnames.yhat))
        c(paste("not", colnames.yhat[1]), colnames.yhat[1])
    else if(!is.null(colnames(fitted.values)))
        c(paste("not", colnames(fitted.values)), colnames(fitted.values))
    else
        last.resort
}
get.class.names <- function(y, yhat, colnames.yhat, fitted.values)
{
    ynames <- paste0("response", seq_len(ncol(yhat)))
    if(length(colnames.yhat) == ncol(yhat))
         class.names <- colnames.yhat
    else if(length(colnames(y)) == ncol(y))
        class.names <- colnames(y)
    else if(length(fitted.values(y)) == fitted.values(y))
        class.names <- colnames(y)
    else
        class.names <- ynames
    # fill in missing names, if necessary
    which. <- which(class.names == "")
    if(length(which.))
        class.names[which.] <- ynames[which.]
    class.names
}
# return names1 but with yhat column names prefixed if necessary
get.prefixed.names <- function(names1, yhat, colnames.yhat, nresponse)
{
    stopifnot(NCOL(yhat) == 1)
    if(!is.null(colnames.yhat) && !is.null(nresponse))
        names1 <- paste(colnames.yhat[nresponse], names1, sep=": ")
    names1
}
# determine the threshold to split classes, a bit of a hack
get.thresh <- function(y, yname)
{
    thresh <- 0
    ymin <- min(y)
    if(ymin == 1)
        thresh <- 1
    if(!is.null(colnames(y)))
        yname <- colnames(y)[1]
    if(ymin == thresh) {
        text.le <- sprintf("%s == %g", yname, thresh)
        text.gt <- sprintf("%s != %g",  yname, thresh)
    } else {
        text.le <- sprintf("%s <= %g", yname, thresh)
        text.gt <- sprintf("%s > %g",  yname, thresh)
    }
    list(thresh=thresh, text.le=text.le, text.gt=text.gt)
}
cannot.plot.this.response <- function(y, yhat, colnames.yhat, type, nlevs, ...)
{
    stop0("cannot plot this kind of response (with predict type=\"",
          type, "\")\n", ...,
          "Additional information:\n  class(observed)=", class(y[1]),
          if(nlevs > 0) paste0(" nlevels(observed)=", nlevs) else "",
          "   ncol(observed)=", NCOL(y),
          if(!is.null(colnames(y)))
              sprintf("   colnames(observed) %s", paste.trunc(colnames(y)))
          else
              "",
          "\n  class(predicted)=", class(yhat[,1])[1],
          "   ncol(predicted)=", NCOL(yhat),
          if(!is.null(colnames.yhat))
              sprintf("   colnames(response) %s", paste.trunc(colnames.yhat))
          else
              "")
}
trace.response.type <- function(trace, type, observed, predicted, ...)
{
    if(trace >= 1) {
        trace2(trace, "\n")
        cat0("observed response: ", observed, "\n",
             "predicted response: ", predicted,
             "   (predict type is \"", type, "\")\n",
             "grouping criterion: ", ..., "\n")
    }
}
draw.legend <- function(densities, degenerate, yhat.per.class, ymax,
                    hist, xlim, col, fill, lty, legend.names,
                    legend.pos, cex.legend, legend.bg, legend.extra)
{
    get.legend.pos <- function()
    {
        # take a stab at positioning the legend correctly --
        # on left or right, away from the highest peak
        pos <- c(0,ymax)
        pos[1] <- xlim[1] # place on left side of graph
        max.left <- 0
        max.right <- 0
        xmid <- xlim[1] + (xlim[2] - xlim[1])/2
        for(iclass in seq_len(nclasses)) {
            if(!degenerate[iclass]) {
                den <- densities[[iclass]]
                if(hist)
                    den$x <- den$x[-1]
                x.left <- (den$x >= xlim[1]) & (den$x <= xmid)
                if(sum(x.left))
                    max.left <- max(max.left, den$y[x.left])
                x.right <- (den$x > xmid) & (den$x <= xlim[2])
                if(sum(x.right))
                    max.right <- max(max.right, den$y[x.right])
            }
        }
        if(max.right < max.left)
            pos[1] <- xlim[1] + (xlim[2] - xlim[1]) / 2.1 # slightly to left of center
        pos
    }
    #--- draw.legend starts here ---

    nclasses <- length(yhat.per.class)
    if(is.null(legend.pos))
        legend.pos <- get.legend.pos()
    else if(length(legend.pos) == 1)
        legend.pos <- c(legend.pos, 0)
    if(length(legend.pos) != 2)
        stop0("length(legend.pos) != 2")
    if(is.null(legend.names))
        legend.names <- names(yhat.per.class)
    if(length(legend.names) < nclasses) {
        warning0("length ", length(legend.names), " of legend.names ",
                 "is less than the number ", nclasses, " of classes")
        legend.names <- repl(legend.names, nclasses)
    }
    else for(iclass in seq_len(nclasses))
        if(degenerate[iclass])
            legend.names[iclass] <- paste(legend.names[iclass], "(not plotted)")
    lwd <- repl(1, nclasses)
    # if the first histogram is filled in, then make its legend lwd bigger
    if(fill[1]==col[1] && fill[1] != "white" && fill[1] != 0)
        lwd[1] <- 4

    if(legend.extra)
        legend.names <- paste0(legend.names, " (", sapply(yhat.per.class, length),
                               " cases)")

    legend(x=legend.pos[1], y=legend.pos[2], legend=legend.names,
           cex=cex.legend, bg=legend.bg, lty=lty, lwd=lwd, col=col)
}
# shade the "error areas" of the density plots

draw.err.col <- function(densities, thresh, col, border, lwd)
{
    den1 <- densities[[1]]
    den2 <- densities[[2]]
    # is reducible error area to the left or to the right?
    # set iden=1 if to the left, iden=2 if to the right
    iden <- den1$y[den1$x >= thresh][1] > den2$y[den2$x >= thresh][1]
    if(is.na(iden)) { # no overlap between classes?
        warning0("no overlap between (first two) classes, ignoring 'err.col' argument")
        return(NULL)
    }
    iden <- if(iden) 2 else 1
    if(length(col) < 2)
        col[2] <- col[1]
    if(length(col) < 3)
        col[3] <- col[iden]
    if(length(border) < 2)
        border[2] <- border[1]
    if(length(border) < 3)
        border[3] <- border[iden]
    if(length(lwd) < 2)
        lwd[2] <- lwd[1]
    if(length(lwd) < 3)
        lwd[3] <- lwd[iden]
    if(is.specified(col[1]) || is.specified(border[1])) {
        # left side of threshold
        matches <- den2$x < thresh
        if(sum(matches)) {
            x <- c(den2$x[matches])
            y <- c(den2$y[matches])
            len <- length(x)
            x[len] <- thresh # close possible tiny gap
            x[len+1] <- thresh
            y[len+1] <- 0
            polygon(x, y, col=col[1], border=border[1], lwd=lwd[1])
        }
    }
    if(is.specified(col[2]) || is.specified(border[2])) {
        # right side of threshold
        matches <- den1$x > thresh
        if(sum(matches)) {
            x <- den1$x[matches]
            y <- den1$y[matches]
            x[1] <- thresh # close possible tiny gap
            len <- length(x)
            x[len+1] <- thresh
            y[len+1] <- 0
            polygon(x, y, col=col[2], border=border[2], lwd=lwd[2])
        }
    }
    if(is.specified(col[3]) || is.specified(border[3])) {
        if(iden == 1) {
            # reducible error, left side of threshold
            # get indices i1 of den1 and i2 of den2 where den1 crosses den2
            i2 <- length(den2$x)
            for(i1 in length(den1$x):1) {
                while(i2 > 1 && den2$x[i2] > den1$x[i1])
                    i2 <- i2 - 1
                if(den1$x[i1] <= thresh && den1$y[i1] >= den2$y[i2])
                    break
            }
            i1 <- den1$x <= thresh & (1:length(den1$x)) >= i1
            i2 <- den2$x <= thresh & (1:length(den2$x)) >= i2
            if(sum(i1) && sum(i2)) {
                # reverse i1 so polygon ends where it starts
                i1 <- rev(((1:length(den1$x))[i1]))
                # close tiny x gap to left of threshhold line
                den1$x[i1][1] <- thresh
                den2$x[i2][sum(i2)] <- thresh
                polygon(x = c(den1$x[i1], den2$x[i2]),
                        y = c(den1$y[i1], den2$y[i2]),
                        col=col[3], border=border[3], lwd=lwd[3])
           }
        } else {
            # reducible error, right side of threshold
            # get indices i1 of den1 and i2 of den2 where den1 crosses den2
            i2 <- 1
            for(i1 in seq_along(den1$x)) {
                while(i2 < length(den2$x) && den2$x[i2] < den1$x[i1])
                    i2 <- i2 + 1
                if(den1$x[i1] >= thresh && den2$y[i2] >= den1$y[i1])
                    break
            }
            i1 <- den1$x >= thresh & (1:length(den1$x)) < i1
            i2 <- den2$x >= thresh & (1:length(den2$x)) < i2
            if(sum(i1) && sum(i2)) {
                # reverse i2 so polygon ends where it starts
                i2 <- rev(((1:length(den2$x))[i2]))
                polygon(x = c(den1$x[i1], den2$x[i2]),
                        y = c(den1$y[i1], den2$y[i2]),
                        col=col[3], border=border[3], lwd=lwd[3])
            }
        }
    }
}
