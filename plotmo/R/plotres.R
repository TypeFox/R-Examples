# plotres.R: plot model residuals

# values for which
W1       <- 1   # model selection
W2CUM    <- 2   # cumulative distribution
W3RESID  <- 3   # residuals vs fitted
W4QQ     <- 4   # qq plot
W5ABS    <- 5   # abs residuals vs fitted
W6SQRT   <- 6   # sqrt abs residuals vs fitted
W7VLOG   <- 7   # abs residuals vs log fitted
W8CUBE   <- 8   # cube root of the squared residuals vs log fitted
W9LOGLOG <- 9   # log abs residuals vs log fitted

# values for vs
V1FITTED   <- 1  # fitted
V2INDEX    <- 2  # obs number
V3RESPONSE <- 3  # response
V4LEVER    <- 4  # leverages

plotres <- function(object = stop("no 'object' argument"),
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

    object.name = quote.deparse(substitute(object)),

    ...)                # passed to predict
{
    init.global.data()
    on.exit(init.global.data()) # release memory on exit
    object # make sure object exists
    trace <- as.numeric(check.integer.scalar(trace, logical.ok=TRUE))
    # Associate the model environment with the object.
    # (This is instead of passing it as an argument to plotmo's data access
    # functions.  It saves a few hundred references to model.env in the code.)
    attr(object, ".Environment") <- get.model.env(object, object.name, trace)
    temp <- plotmo_prolog(object, object.name, trace, ...)
        object  <- temp$object
        my.call <- temp$my.call
    if(!is.numeric(which) || !is.vector(which) || anyNA(which) ||
            any(which != floor(which)) || any(which < 1) || any(which > W9LOGLOG))
        which.err()
    info        <- check.boolean(info)
    standardize <- check.boolean(standardize)
    delever     <- check.boolean(delever)
    level       <- check.level.arg(level, zero.ok=TRUE)
    smooth.col  <- get.smooth.col(smooth.col, ...)
    grid.col    <- dot("col.grid", DEF=grid.col, ...)
    if(is.specified(grid.col) && is.logical(grid.col) && grid.col) # grid.col=TRUE
        grid.col <- "lightgray"
    check.integer.scalar(nresponse, min=1, na.ok=TRUE, logical.ok=FALSE, char.ok=TRUE)
    # set random seed for reproducibility if jitter or isubset is used
    rnorm(1) # seems to be necessary to make .Random.seed available
    old.seed <- .Random.seed
    set.seed(2015)
    on.exit(set.seed(old.seed), add=TRUE)
    temp <- get.plotres.data(object, object.name, which, standardize, delever,
                             level, versus, id.n, labels.id,
                             trace, npoints, type, nresponse, ...)
    nresponse <- temp$nresponse   # col index in the response (converted from NA if necessary)
    resp.name <- temp$resp.name   # used only in automatic caption, may be NULL
    type      <- temp$type        # always a string (converted from NULL if necessary)
    rinfo     <- temp$rinfo       # resids, scale, name, etc.
    vinfo     <- temp$vinfo       # versus.mat, icolumns, nversus, etc.
    fitted    <- temp$fitted      # n x 1 numeric matrix, colname is "Fitted"
    which     <- temp$vinfo$which # plots we don't want will have been removed
    id.n      <- temp$id.n        # forced to zero if row indexing changed
    npoints   <- temp$npoints     # special values have been converted
    rsq       <- temp$rsq         # r-squared on the training data

    possible.biglm.warning(object, trace)

    nfigs <- length(which) * length(vinfo$icolumns)
    if(nfigs == 0) {
        if(trace >= 0)
            warning0("plotres: nothing to plot")
        return(invisible())
    }
    do.par <- check.do.par(do.par, nfigs) # do.par is 0, 1, or 2
    # Prepare caption --- we need it now for do.par() but
    # can only display it later after at least one plot.
    caption <- get.caption(nfigs, do.par, caption, resp.name, type,
                           object$call, object.name, my.call)
    if(do.par) {
        oldpar <- par(no.readonly=TRUE)
        do.par(nfigs = nfigs, caption=caption, main1=dot("main", ...),
               xlab1 = dot("xlab", DEF=NULL, ...), # used only for margin spacing
               ylab1 = dot("ylab", DEF=NULL, ...), # ditto
               trace = trace,
               def.font.main = 1, # for compat with lm.plot
               nlines.in.main = # nbr of lines in main is needed for margins
                nlines.in.main(object=object, which=which, versus=versus,
                    standardize=standardize, delever=delever, level=level, ...),
               ...)
        if(do.par == 1)
            on.exit(par(oldpar), add=TRUE)
    } else { # do.par=FALSE
        oldpar <- do.par.dots(..., trace=trace)
        if(length(oldpar))
            on.exit(do.call(par, oldpar), add=TRUE)
    }
    # TODO if length(which) == 1 then pass xlim and ylim to the specified plot
    #      else if(1 %in% which) pass xlim and ylim to the which=1 plot only
    #      else pass xlim and ylim to the residuals plots

    force.auto.resids.xylim <- (W1 %in% which)

    if(any(which == W1)) {
        plotted <- plot_w1(object=object, which=which, info=info,
            standardize=standardize, delever=delever, level=level, versus=versus,
            id.n=id.n, labels.id=rinfo$labs, smooth.col=smooth.col,
            grid.col=grid.col, do.par=do.par,
            # must do caption here if will not call plot1 later
            caption=if(all(which == W1)) caption else "", trace=trace,
            npoints=npoints, center=center, type=type, nresponse=nresponse,
            object.name=object.name,
            ...)
        which <- which[which != W1]
        if(length(which) == 0 && !plotted && trace >= 0)
            warning0("plotres: nothing to plot")
    }
    if(length(which) == 0)
        return(invisible())

    # we do this after the w1 call so we pass NULL to w1 if labels.id were NULL
    if(is.null(rinfo$labs))
        rinfo$labs <- paste(1:length(rinfo$resids))

    # We plot only the residuals in iresids, but use all the
    # residuals for calculating densities (where "all" actually means
    # a maximum of NMAX cases, see the previous call to get.isubset).
    #
    # The "use.all=(nversus == V4LEVER)" keeps things easier later
    # for leverage plots, but it would be nice to not have to use it.

    iresids <- get.isubset(rinfo$resids, npoints, id.n,
                           use.all=(vinfo$nversus == V4LEVER), rinfo$scale)

    xlim <- dot("xlim", DEF=NULL, ...)

    for(icolumn in vinfo$icolumns) {
        for(iwhich in seq_along(which)) {
            if(which[iwhich] == W2CUM)
                plotmo_cum(rinfo=rinfo, info=info, nfigs=nfigs, add=FALSE,
                           cum.col1=NA, grid.col=grid.col, jitter=0, ...)
            else if(which[iwhich] == W4QQ)
                plotmo_qq(rinfo=rinfo, info=info, nfigs=nfigs,
                          grid.col=grid.col, smooth.col=smooth.col,
                          id.n=id.n, iresids=iresids, npoints=npoints, ...)
            else
                plotresids(object=object, which=which[iwhich],
                    info=info, standardize=standardize, level=level,
                    # versus1 is what we plot along the x axis, a vector
                    versus1=vinfo$versus.mat[, icolumn],
                    id.n=id.n, smooth.col=smooth.col,
                    grid.col=grid.col, jitter=jitter,
                    npoints=npoints, center=center,
                    type=type,
                    fitted=fitted, rinfo=rinfo, rsq=rsq, iresids=iresids,
                    nversus=vinfo$nversus,
                    colname.versus1=colnames(vinfo$versus.mat)[icolumn],
                    force.auto.resids.xylim=force.auto.resids.xylim,
                    ...)
        }
    }
    draw.caption(caption, ...)
    if(trace >= 1)
        printf("\ntraining rsq %.2f\n", rsq)
    invisible()
}
which.err <- function()
{
    stop0("Bad value for which\n",
          "Allowed values for which:\n",
          "  1  Model\n",
          "  2  Cumulative distribution\n",
          "  3  Residuals vs fitted\n",
          "  4  QQ plot\n",
          "  5  Abs residuals vs fitted\n",
          "  6  Sqrt abs residuals vs fitted\n",
          "  7  Abs residuals vs log fitted\n",
          "  8  Cube root of the squared residuals vs log fitted\n",
          "  9  Log abs residuals vs log fitted")
}
versus.err <- function()
{
    stop0("versus must be an integer or a string:\n",
          "  1    fitted (default)\n",
          "  2    observation numbers\n",
          "  3    response\n",
          "  4    leverages\n",
          "  \"\"   predictors\n",
          "  \"b:\" basis functions")
}
nlines.in.main <- function(object, which, versus, standardize, delever, level, ...)
{
    main <- dot("main", ...)

    if(is.specified(main)) # user specified main?
        return(nlines(main))

    # auto main
    has.extra.line <- # conservative guess if main will have two lines
        standardize || delever || level ||
        any(which %in% W6SQRT:W9LOGLOG) ||
        (versus %in% V4LEVER) || is.character(versus)

    max(1 + has.extra.line, nlines.in.w1.main(object))
}
