# w1.R: plotres functions for the which=1 plot

plot_w1 <- function(object,
    which, # currently used only to get the total nbr of plots (for xlim and ylim)
    info, standardize, delever, level, versus,
    id.n, labels.id, smooth.col, grid.col,
    do.par, caption, trace,
    npoints, center,
    type, nresponse,
    object.name,
    SHOWCALL=NA, # this is here to absorb SHOWCALL from dots
    ...)
{
    call.w1 <- function(FUNC, ...)
    {
        fname <- trunc.deparse(substitute(FUNC))
        call.dots(FUNC=FUNC, PREFIX="w1.", DROP="*", KEEP="PREFIX",
                  TRACE=trace >= 1,
                  FNAME=fname, force.anon=object, ...)
    }
    plotted = TRUE
    if(inherits(object, "earth"))
        call.earth.modsel(object=object, trace=trace, grid.col=grid.col, ...)
    else if(inherits(object, "mars")) { # mda::mars, convert first to an earth model
        if(trace)
            printf("calling mars.to.earth (needed for the model selection plot)\n")
        earth.mod <- earth::mars.to.earth(object, trace=trace >= 2)
        earth.mod <- update(earth.mod, trace=trace >= 2)
        call.earth.modsel(object=earth.mod, trace=trace, grid.col=grid.col, ...)
    } else if(inherits(object, "lm")) {
        # check that the model supports hatvalues(), needed for versus=V4LEVER.
        if(is.try.err(try(hatvalues(object), silent=TRUE)))
            plotted <- FALSE
        else {
            # do a recursive call to plotres to plot the residuals versus leverage plot
            if(trace >= 1)
                printf(
"plotres(object, which=3, versus=4, ...)  (recursive call for leverage plot)\n")
            plotres(object=object, which=W3RESID, info=info,
                standardize=standardize, delever=delever, level=level,
                versus=V4LEVER,
                id.n=id.n, labels.id=labels.id, smooth.col=smooth.col,
                grid.col=grid.col,
                do.par=FALSE, caption=caption,
                trace=if(trace==1) 0 else trace,
                npoints=npoints, center=center,
                type=type, nresponse=nresponse,
                object.name=object.name,
                ...)
        }
    } else if(inherits(object, "rpart")) {
        if(requireNamespace("rpart.plot", quietly=TRUE))
            call.w1(rpart.plot::rpart.plot, ...)
        else {
            printf("Use the \"rpart.plot\" package for better rpart plots\n")
            plot(object, compress=TRUE, uniform=TRUE)
            text(object, xpd=NA)
        }
    } else if(inherits(object, "tree")) {
        call.w1(graphics::plot, def.type="uniform", ...)
        n <- nrow(object$frame)
        def.cex <- if(n < 8) 1 else if(n < 20) .9 else .8
        call.w1(graphics::text, def.pretty=3, def.digits=3, def.cex=def.cex, ...)
    } else if(inherits(object, "randomForest"))
        call.w1(graphics::plot, ...,
                def.main=dot("main", DEF="Error vs Number of Trees", ...))
    else if(inherits(object, "gbm")) {
        # don't allow w1.n.trees argument, except w1.n.trees=NA
        predict.n.trees <- dot("predict.n.trees", DEF=object$n.trees,  ...)
        w1.n.trees      <- dot("w1.n.trees",      DEF=predict.n.trees, ...)
        if(!is.na(w1.n.trees) && w1.n.trees != predict.n.trees) {
            if(is.na(dot("predict.n.trees", EX=0, ...)))
                stop0("w1.n.trees is not allowed (please use predict.n.trees)")
            else
                stop0("w1.n.trees is not allowed")
        }
        call.w1(plot.gbmx, w1.n.trees=w1.n.trees, ...)
    } else if(inherits(object, "cosso"))
        call.w1(graphics::plot, def.M=2, ...)
    else if(inherits(object, c("glmnet", "multnet", "mrelnet")))
        call.w1(plot.glmnetx, def.xvar="rlambda", def.grid.col=grid.col,
                force.s=attr(object, "s"), force.nresponse=nresponse, ...)
    else if(inherits(object, c("lars", "sparsenet", "cv.glmnet")))
        call.w1(graphics::plot, ...)

    # # TODO commented out because plot.nn uses grid graphics
    # #      which doesn't coexist with base graphics
    # } else if(inherits(object, "nn")) { # neuralnet package
    #     rep <- dot("w1.rep", DEF="best", ...)
    #     if(is.null(rep))
    #         stop0("rep=NULL is not allowed here for plot.nn ",
    #               "(because it invokes dev.new)")
    #     call.w1(plot.nn, def.rep=rep, ...)

    else
        plotted <- FALSE

   draw.caption(caption, ...)
   plotted
}
# Note that by specifying col and lty in the arg list we drop
# them from dots passed to earth_plotmodsel, else get
# 'col' matches both the 'col.rsq' and 'col.grsq' arguments.
# TODO call.dot should be able to do this dropping for us but currently can't
call.earth.modsel <- function(object, trace, grid.col, col=NA, lty=NA, ...)
{
    call.dots(earth::earth_plotmodsel, PREFIX="w1.",
             DROP="*", KEEP="PREFIX,PLOT.ARGS,PLOTMO.ARGS",
             trace=trace >= 1,
             force.x=object, grid.col=grid.col, ...)
}
nlines.in.w1.main <- function(object)
{
    1
}
