
### Some obsolete experimental stuff. Remove at some point.


### Now, the question is, (1) whether to have the method for data
### frames (which is what's of immediate interest) here or in lme4.
### It doesn't really make sense to have it in lme4, since we won't be
### defining a new class.  But most of the features would be borrowed
### from nlme (2) whether to have methods for any other class (like
### 'lm' maybe)



### data frame method: gplotArgs.data.frame(x, <special args>, ...)

## 4 levels of information:
##    o attr(x, "ginfo")
##    o <special args>
##    o attr(x, "gplot.args")
##    o ...
## the first 2 determine a reasonable list, overridden by the 3rd and
## then the 4th



## should we have a non-generic "gplotArgs<-" for setting the
## "gplot.args" attribute?



## First, we need some good default panel functions.  These need to
## depend on

## (1) the type of display formula -
##      o factor ~ numeric   [ model: numeric ~ 1       | factor ]
##      o numeric ~ numeric  [ model: numeric ~ numeric | factor ]
##      o  ~ numeric         [ model:  factor ~ numeric | factor ] NEW


panel.df.fn <-               # factor ~ numeric
    function(x, y, groups = NULL, ...)
{
    panel.dotplot(x=x, y=y, groups = groups, ...)
}

panel.df.nn <-               # numeric ~ numeric
    function(x, y, lines = all(type == "p"), type = "p", ...)
{
    panel.xyplot(x, y, type = type, ...)
    if (lines)
    {
        y.avg <- tapply(y, x, mean) # lines through average y
        y.avg <- y.avg[!is.na(y.avg)]
        if (length(y.avg) > 0)
        {
            xvals <- as.numeric(names(y.avg))
            ord <- order(xvals)
            panel.xyplot(xvals[ord], y.avg[ord], type = "l", ...)
        }
    }
}



panel.df <- # combines above 2, should be called
    function(x, y = NULL, groups = NULL, grid = TRUE, ...)
{
    yNull <- is.null(y)
    groupsNull <- is.null(groups)
    xFactor <- is.factor(x)
    yFactor <- is.factor(y)

    if (yNull) ## for factor response, may not be in use yet
    {
        if (groupsNull) panel.densityplot(x = x, ...)
        else panel.superpose(x = x, groups = groups,
                             panel.groups = "panel.densityplot", ...)
    }
    else if (!xFactor && !yFactor) ## numeric ~ numeric
    {
        if (grid) panel.grid()
        if (groupsNull) panel.df.nn(x = x, y = y, ...)
        else panel.superpose(x = x, y = y, groups = groups,
                             panel.groups = panel.df.nn, ...)
    }
    else if (yFactor && !xFactor)
    {
        panel.df.fn(x = x, y = y, groups = groups, ...)
    }
    else stop("can't handle both x and y being factors yet")
    invisible()
}


## display formula

## generally, the display is determined by the arguments 'formula' and
## 'groups'.  In the nlme scheme, the formula associated with the data
## is a model formula, not a display formula.  The display is further
## determined by the 'outer' and 'inner' arguments.  

## Let's try to outline a strategy for this.

## attr(x, "ginfo") can contain things used in nlme's groupedData.

## formula can be of the form 'y ~ x | id' (only one variable in 'id',
## although it can be an interaction). inner=~a becomes default
## grouping factor. outer=~b+c used for ordering levels of id by
## values of y.  outer can be also an argument to
## gplotArgs.data.frame.  If TRUE, it's equivalent to outer =
## ginfo$outer, or it could be =~e+f.  In either case, it becomes the
## conditioning variables, id becomes groups (inner would then be
## ignored)

## In the formula itself, if x is a factor, the display formula
## becomes id ~ y, with groups = x (unless x = 1, then no groups)



## .typeInDF <- function(x, data)
## {
##     if (x == "1") "one"
##     else if (is.factor(data[[x]])) "factor"
##     else if (is.numeric(data[[x]])) "numeric"
##     else stop(paste("don't recognize", class(data[[x]])))
## }



## should work with expressions like log(height) as well
.typeInDF <- function(x, data)
{
    if (all(is.na(x))) return(NA)
    x <- eval(parse(text = x), data, parent.frame())
    if (length(x) == 1 && x == 1) "one"
    else if (is.factor(x)) "factor"
    else if (is.numeric(x)) "numeric"
    else stop(paste("don't recognize", class(x)))
}




## S3 method: underlies groupedData() in lme4

## FIXME: the following doesn't work in the docs.  What's the recommended way?
## \method{"gplotArgs<-"}{data.frame}(x, value)



"gplotArgs<-.data.frame" <- 
    function(x, value)
{
    if (!is.list(value)) stop("assigned value must be list")
    process.args <-
        function(formula, 
                 order.groups = TRUE,
                 FUN = function(x, ...) max(x, na.rm = TRUE),
                 outer = NULL, inner = NULL,
                 labels = list(), units = list(), ...)
        {
            list(ginfo =
                 list(formula = formula,
                      order.groups = order.groups,
                      FUN = FUN,
                      outer = outer,
                      inner = inner,
                      labels = labels,
                      units = units),
                 dots = list(...))
        }
    pargs <- do.call(process.args, value)
    attr(x, "ginfo") <- pargs$ginfo
    if (length(pargs$dots) > 0)
        attr(x, "gplot.args") <- pargs$dots
    x
}


## S3 method

## basic idea: get defaults from "ginfo" attribute, then overwrite by
## "gplot.args" attribute, followed by ...

gplotArgs.data.frame <-
    function(x, display.formula, outer = FALSE, inner = FALSE,
             groups = NULL,
             ...,
             subset = TRUE)
{

    ## The final result should contain only standard trellis args,
    ## with a special component plotFun, and optionally a
    ## display.formula, which overrides formula.  However, for data
    ## frames, there are some other issues.

    ## The "ginfo" attribute can only contain info traditionally in
    ## nlme groupedData objects.  This, along with explicit arguments
    ## here, will be used to create a list of trellis args.  These can
    ## be overridden by the "gplot.args" attribute, and then by
    ## ... here.


    ginfo <- attr(x, "ginfo")
    gplot.args <- attr(x, "gplot.args")
    ## equivalent to default method if ginfo is NULL
    if (is.null(ginfo) || !is.list(ginfo))
        return(updateList(gplot.args, list(...)))

    ## First (longish) task: determine display formula and groups

    ## groups is by far the most irritating thing to handle.  The only
    ## options I can think of: (1) evaluate groups here and pass it on
    ## and (2) change lattice to allow groups to be a formula.  Use
    ## (1) for now

    groups <- eval(substitute(groups), x, parent.frame()) # typically NULL

    ## subset poses a similar problem

    subset <- eval(substitute(subset), x, parent.frame())


    if (missing(display.formula))
        display.formula <- gplot.args$display.formula

    ## major step: get the display formula, but only if it's NULL
    model.formula <- ginfo$formula
    if (!is.null(display.formula))
        ## no point in jumping through hoops
    {
        ans <- list(display.formula = display.formula)
    }
    else if (!is.null(model.formula)) ## the interesting stuff
    {
        vars <-
            list(resp = .responseName(model.formula),
                 cov = .covariateName(model.formula),
                 grp = .groupsName(model.formula))




        if (ginfo$order.groups)
        {
            ## reorder grp based on values of resp

            if (is.null(ginfo$FUN)) ginfo$FUN <- function(x, ...) max(x, na.rm = TRUE)

            respVar <- vars$resp
            grpVar <- vars$grp

            scores <- tapply(x[[respVar]], x[[grpVar]], ginfo$FUN)

            if (inherits(ginfo$outer, "formula"))
            {
                outerVar <- .covariateName(ginfo$outer)
                outer.unique <- tapply(x[[outerVar]], x[[grpVar]], unique)
                ord <- order(outer.unique, scores)
            }
            else
                ord <- order(scores)

            x[[grpVar]] <- factor(x[[grpVar]], levels = names(scores)[ord])
        }






        ## display formula may be further modified by inner and outer.
        ## How does that affect rest of the calculations?

        if (is.logical(outer) && outer) outer <- ginfo$outer
        if (is.logical(inner) && inner) inner <- ginfo$inner

        ## both cannot happen. outer makes outer the conditioning
        ## variables, normal grp becomes groups. inner behaves as
        ## groups. outer takes precedence.

        if (inherits(outer, "formula"))
        {
            if (is.null(groups)) groups <- as.formula(paste("~", vars$grp))
            vars$grp <- .covariateName(outer)
        }
        else if (inherits(inner, "formula"))
        {
            ## FIXME: this may not be the right thing to do
            if (is.null(groups)) groups <- inner
        }
        


        varTypes <-
            lapply(vars, .typeInDF, data = x)

        ## Next step depends on varTypes

        ## case 1: cov = "one" - grp ~ resp
        ## case 2: cov = "factor" - grp ~ resp, groups = cov
        ## case 1: cov = "numeric" - resp ~ cov | grp

        fc <- switch(varTypes$cov,
                     one = paste(vars$grp, "~", vars$resp),
                     factor = paste(vars$grp, "~", vars$resp),
                     numeric =
                     paste(vars$resp, "~", vars$cov, "|", vars$grp))
        ans <-
            list(display.formula = as.formula(fc))

        if (varTypes$cov == "factor" && is.null(groups))
            groups <-
                eval(parse(text = vars$cov), x, parent.frame())
    }
    else 
    {
        stop("no formula!")
    }





    ## determine default plot function based on display.formula

    dvars <-
        list(resp = .responseName(ans$display.formula),
             cov = .covariateName(ans$display.formula),
             grp = .groupsName(ans$display.formula)) ## NA is none
    dvarTypes <-
        lapply(dvars, .typeInDF, data = x)
    if (dvarTypes$resp == "numeric" && dvarTypes$cov == "numeric")
        plotFun.constructed <- "xyplot"
    else if (dvarTypes$resp == "factor" && dvarTypes$cov == "numeric")
        plotFun.constructed <- "dotplot"
    else {
        ## str(dvarTypes)
        stop("unsupported combination")
    }

    ## other stuff in ginfo?

    ylab.constructed <-
        if ("labels" %in% names(ginfo) && dvars$resp %in% names(ginfo$labels))
            ginfo$labels[[dvars$resp]]
        else dvars$resp
    xlab.constructed <-
        if ("labels" %in% names(ginfo) && dvars$cov %in% names(ginfo$labels))
            ginfo$labels[[dvars$cov]]
        else dvars$cov
    if ("units" %in% names(ginfo) && dvars$resp %in% names(ginfo$units))
        ylab.constructed <- paste(ylab.constructed, ginfo$units[[dvars$resp]])
    if ("units" %in% names(ginfo) && dvars$cov %in% names(ginfo$units))
        xlab.constructed <- paste(xlab.constructed, ginfo$units[[dvars$cov]])


    ans <-
        updateList(ans,
                   list(plotFun = plotFun.constructed,
                        data = x,
                        panel = panel.df,
                        groups = groups,
                        subset = subset,
                        xlab = xlab.constructed, ylab = ylab.constructed,
                        aspect = if (plotFun.constructed == "xyplot") "xy" else "fill",
                        auto.key =
                        switch(plotFun.constructed,
                               xyplot = list(points = FALSE, lines = TRUE, space = "right"),
                               dotplot = list(points = TRUE, space = "right"))))


    
    if (!is.null(gplot.args))
        ans <- updateList(ans, gplot.args) ## leave this out?
    updateList(ans, list(...))
}




