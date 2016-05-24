

### Copyright (C) 2005 Deepayan Sarkar
### <Deepayan.Sarkar@R-project.org>, Douglas Bates
### <Douglas.Bates@R-project.org>.  See file COPYING for license
### terms.





### unexported helper function to obtain a valid subset argument.
### Mostly to ensure that we don't up with a huge vector of integer
### indices in the trivial cases.



thinned.indices <- function(object, n = NROW(object), start = 1, thin = 1)
{
    if (is.mcmc(object) &&
        (start * thin != 1) &&
        !all(mcpar(object)[-2] == 1))
        warning("mcmc object is already thinned")
    if (start < 1) stop("Invalid start")
    else if (start == 1)
    {
        if (thin < 1) stop("Invalid thin")
        else if (thin == 1) TRUE
        else rep(c(TRUE, FALSE), c(1, thin-1))
    }
    else
    {
        if (thin < 1) stop("Invalid thin")
        else if (thin == 1) -seq(length = start-1)
        else start + thin * (0:(floor(n - start) / thin))
    }
}



## most functions will have methods for mcmc and mcmclist objects.  In
## the second case, another grouping variable is added.  By default,
## this will be used for ``grouped displays'' (when possible), but
## could also be used for conditioning by setting groups=FALSE





## levelplot, analog of plot.crosscorr.  Defaults changed to match
## those of plot.crosscorr as much as possible.


levelplot.mcmc <-
    function(x, data = NULL,
             main = attr(x, "title"),
             start = 1, thin = 1,
             ...,
             xlab = "", ylab = "",
             cuts = 10,
             at = do.breaks(c(-1.001, 1.001), cuts),
             col.regions = topo.colors(100),
             subset = thinned.indices(x, start = start, thin = thin))
{
    if (!is.R()) {
      stop("This function is not yet available in S-PLUS")
    }
    cormat <- cor(x[subset, ])
    cormat <- cormat[, rev(seq(length = ncol(cormat)))]
    levelplot(cormat,
              main = main, ...,
              cuts = cuts, at = at,
              col.regions = col.regions,
              xlab = xlab, ylab = ylab)
}


## the mcmc.list method wouldn't do any grouping (so should have
## outer=TRUE by default).  It hasn't been written yet because 




## The splom (FIXME: not yet written) method may be more useful

## in progress, unexported.  Planning to make it be like levelplot on
## the lower diagonal (maybe ellipses instead of plain boxes) and
## normal splom on the upper diagonal.  Not much point in having tick
## marks.  Names in the middle save space, unlike in levelplot which
## stupidly shows correlation=1 on the diagonal.


splom.mcmc <-
    function(x, data = NULL,
             main = attr(x, "title"),
             start = 1, thin = 1,
             as.matrix = TRUE,
             xlab = "", ylab = "",
             cuts = 10,
             at = do.breaks(c(-1.001, 1.001), cuts),
             col.regions = topo.colors(100),
             ...,
             pscales = 0,
             subset = thinned.indices(x, start = start, thin = thin))
{
##     cormat <- cor(x[subset, ])
##     cormat <- cormat[, rev(seq(length = ncol(cormat)))]
    if (!is.R()) {
      stop("This function is not yet available in S-PLUS")
    }
    splom(as.data.frame(x[subset, ]),
          as.matrix = as.matrix,
          main = main, ...,
          pscales = pscales,
          cuts = cuts, at = at,
          lower.panel = function(x, y, ...) {
              corval <- cor(x, y)
              grid::grid.text(lab = round(corval, 2))
          },
          col.regions = col.regions,
          xlab = xlab, ylab = ylab)
}








### methods for densityplot (mcmc and mcmc.list)



densityplot.mcmc <-
    function(x, data = NULL,
             outer, aspect = "xy",
             default.scales = list(relation = "free"),
             start = 1, thin = 1,
             main = attr(x, "title"),
             xlab = "",
             plot.points = "rug",
             ...,
             subset = thinned.indices(x, start = start, thin = thin))
{
    if (!is.R()) {
      stop("This function is not yet available in S-PLUS")
    }
    if (!missing(outer)) warning("specification of outer ignored")
    data <- as.data.frame(x)
    form <-
        as.formula(paste("~",
                         paste(lapply(names(data), as.name),
                               collapse = "+")))

    
### The following is one possible approach, but it does not generalize
### to mcmclist objects.


##     densityplot(form, data = data,
##                 outer = outer,
##                 aspect = aspect,
##                 default.scales = default.scales,
##                 main = main,
##                 xlab = xlab,
##                 plot.points = plot.points,
##                 subset = eval(subset),
##                 ...)

    
### This one does, with the only downside I can think of being that
### subscripts, if used, will give indices in subsetted data, not
### original.  But that's true even if the original mcmc object was
### itself already thinned.

    densityplot(form, data = data[subset, , drop=FALSE],
                outer = TRUE,
                aspect = aspect,
                default.scales = default.scales,
                main = main,
                xlab = xlab,
                plot.points = plot.points,
                ...)
}


densityplot.mcmc.list <-
    function(x, data = NULL,
             outer = FALSE, groups = !outer,
             aspect = "xy",
             default.scales = list(relation = "free"),
             start = 1, thin = 1,
             main = attr(x, "title"),
             xlab = "",
             plot.points = "rug",
             ...,
             subset = thinned.indices(x[[1]], start = start, thin = thin))
{
    if (!is.R()) {
      stop("This function is not yet available in S-PLUS")
    }
    if (groups && outer) warning("'groups=TRUE' ignored when 'outer=TRUE'")
    datalist <- lapply(x, function(x) as.data.frame(x)[subset, ,drop=FALSE])
    data <- do.call("rbind", datalist)
    form <-
        if (outer)
            as.formula(paste("~",
                             paste(lapply(names(data), as.name),
                                   collapse = "+"),
                             "| .run"))
##             as.formula(paste("~",
##                              paste(names(data),
##                                    collapse = "+"),
##                              "| .run"))
        else
            as.formula(paste("~",
                             paste(lapply(names(data), as.name),
                                   collapse = "+")))
##             as.formula(paste("~",
##                              paste(names(data),
##                                    collapse = "+")))
    ##data[["index"]] <- seq(length = nrow(x[[1]]))[subset]
    .run <- gl(length(datalist), nrow(datalist[[1]]))
    if (groups && !outer)
        densityplot(form, data = data,
                    outer = TRUE,
                    groups = .run,
                    aspect = aspect,
                    default.scales = default.scales,
                    main = main,
                    xlab = xlab,
                    plot.points = plot.points,
                    ...)
    else
        densityplot(form, data = data,
                    outer = TRUE,
                    aspect = aspect,
                    default.scales = default.scales,
                    main = main,
                    xlab = xlab,
                    plot.points = plot.points,
                    ...)
}


### methods for qqmath (mcmc and mcmc.list)




qqmath.mcmc <-
    function(x, data = NULL,
             outer, aspect = "xy",
             default.scales = list(y = list(relation = "free")),
             prepanel = prepanel.qqmathline,
             start = 1, thin = 1,
             main = attr(x, "title"),
             ylab = "",
             ...,
             subset = thinned.indices(x, start = start, thin = thin))
{
    if (!is.R()) {
      stop("This function is not yet available in S-PLUS")
    }
    if (!missing(outer)) warning("specification of outer ignored")
    data <- as.data.frame(x)
    form <-
        as.formula(paste("~",
                         paste(lapply(names(data), as.name),
                               collapse = "+")))
    qqmath(form, data = data[subset, ,drop=FALSE],
           outer = TRUE,
           aspect = aspect,
           prepanel = prepanel,
           default.scales = default.scales,
           main = main,
           ylab = ylab,
           ...)
}


qqmath.mcmc.list <-
    function(x, data = NULL,
             outer = FALSE, groups = !outer,
             aspect = "xy",
             default.scales = list(y = list(relation = "free")),
             prepanel = prepanel.qqmathline,
             start = 1, thin = 1,
             main = attr(x, "title"),
             ylab = "",
             ...,
             subset = thinned.indices(x[[1]], start = start, thin = thin))
{
    if (!is.R()) {
      stop("This function is not yet available in S-PLUS")
    }
    if (groups && outer) warning("'groups=TRUE' ignored when 'outer=TRUE'")
    datalist <- lapply(x, function(x) as.data.frame(x)[subset, , drop=FALSE])
    data <- do.call("rbind", datalist)
    form <-
        if (outer)
            as.formula(paste("~",
                             paste(lapply(names(data), as.name),
                                   collapse = "+"),
                             "| .run"))
        else
            as.formula(paste("~",
                             paste(lapply(names(data), as.name),
                                   collapse = "+")))
    ##data[["index"]] <- seq(length = nrow(x[[1]]))[subset]
    ##data[[".run"]] <- gl(length(datalist), nrow(datalist[[1]]))
    .run <- gl(length(datalist), nrow(datalist[[1]]))
    if (groups && !outer)
        qqmath(form, data = data,
               outer = TRUE,
               groups = .run,
               aspect = aspect,
               prepanel = prepanel,
               default.scales = default.scales,
               main = main,
               ylab = ylab,
               ...)
    else
        qqmath(form, data = data,
               outer = TRUE,
               aspect = aspect,
               default.scales = default.scales,
               main = main,
               ylab = ylab,
               ...)
}



### methods for xyplot (mcmc and mcmc.list)



xyplot.mcmc <-
    function(x, data = NULL,
             outer, layout = c(1, nvar(x)),
             default.scales = list(y = list(relation = "free")),
             type = 'l',
             start = 1, thin = 1,
             xlab = "Iteration number",
             ylab = "", 
             main = attr(x, "title"),
             ...,
             subset = thinned.indices(x, start = start, thin = thin))
{
    if (!is.R()) {
      stop("This function is not yet available in S-PLUS")
    }
    if (!missing(outer)) warning("specification of outer ignored")
    data <- as.data.frame(x)
    form <- eval(parse(text = paste(paste(lapply(names(data), as.name),
                       collapse = "+"), "~.index")))
    data[[".index"]] <- time(x)
    xyplot(form, data = data[subset, ],
           outer = TRUE,
           layout = layout,
           default.scales = default.scales,
           type = type,
           xlab = xlab,
           ylab = ylab,
           main = main,
           ...)
}



xyplot.mcmc.list <-
    function(x, data = NULL,
             outer = FALSE, groups = !outer,
             aspect = "xy", layout = c(1, nvar(x)),
             default.scales = list(y = list(relation = "free")),
             type = 'l',
             start = 1, thin = 1,
             xlab = "Iteration number",
             ylab = "",
             main = attr(x, "title"),
             ...,
             subset = thinned.indices(x[[1]], start = start, thin = thin))
{
    if (!is.R()) {
      stop("This function is not yet available in S-PLUS")
    }
    if (groups && outer) warning("'groups=TRUE' ignored when 'outer=TRUE'")
    datalist <- lapply(x, function(x) as.data.frame(x)[subset,,drop=FALSE])
    data <- do.call("rbind", datalist)
    form <-
        if (outer)
            eval(parse(text = paste(paste(lapply(names(data), as.name),
                       collapse = "+"), "~.index | .run")))
        else
            eval(parse(text = paste(paste(lapply(names(data), as.name),
                       collapse = "+"), "~.index")))
##     form <-
##         if (outer)
##             as.formula(paste(paste(names(data),
##                                    collapse = "+"),
##                              "~ index | .run"))
##         else
##             as.formula(paste(paste(names(data),
##                                    collapse = "+"),
##                              "~ index"))
    data[[".index"]] <- time(x)
    .run <- gl(length(datalist), nrow(datalist[[1]]))
    if (groups && !outer)
        xyplot(form, data = data,
               outer = TRUE,
               layout = layout,
               groups = .run,
               default.scales = default.scales,
               type = type,
               main = main,
               xlab = xlab,
               ylab = ylab,
               ...)
    else
        xyplot(form, data = data,
               outer = TRUE,
               layout = layout,
               default.scales = default.scales,
               type = type,
               main = main,
               xlab = xlab,
               ylab = ylab,
               ...)
}








### methods for acfplot (mcmc and mcmc.list)



panel.acfplot <-
    function(..., groups = NULL)
{
    if (!is.R()) {
      stop("This function is not yet available in S-PLUS")
    }
    reference.line <- trellis.par.get("reference.line")
    panel.abline(h = 0,
                 col = reference.line$col, 
                 lty = reference.line$lty, 
                 lwd = reference.line$lwd, 
                 alpha = reference.line$alpha)
    if (is.null(groups)) panel.xyplot(...)
    else panel.superpose(..., groups = groups)
}




acfplot <- function(x, data, ...)
    UseMethod("acfplot")



acfplot.mcmc <-
    function(x, data = NULL,
             outer,
             prepanel = function(x, y, ...) list(ylim= c(-1, 1) * max(abs(y[-1]))),
             panel = panel.acfplot,
             type = "h",
             aspect = "xy",
             start = 1, thin = 1,
             lag.max = NULL,
             ylab = "Autocorrelation",
             xlab = "Lag",
             main = attr(x, "title"),
             ...,
             subset = thinned.indices(x, start = start, thin = thin))
{
    if (!is.R()) {
      stop("This function is not yet available in S-PLUS")
    }
    if (!missing(outer)) warning("specification of outer ignored")
    getAcf <- function(x, lag.max)
    {
        as.vector(acf(x, lag.max = lag.max, plot = FALSE)$acf)
    }
    data <- as.data.frame(apply(as.matrix(x)[subset, ,drop=FALSE], 2, getAcf, lag.max = lag.max))
    form <-
        eval(parse(text = paste(paste(lapply(names(data), as.name),
                   collapse = "+"), "~.lag")))
    data[[".lag"]] <- seq(length = nrow(data))
    xyplot(form, data = data,
           outer = TRUE,
           prepanel = prepanel,
           panel = panel,
           type = type,
           aspect = aspect,
           xlab = xlab,
           ylab = ylab,
           main = main,
           ...)
}




acfplot.mcmc.list <-
    function(x, data = NULL,
             outer = FALSE, groups = !outer,
             prepanel = function(x, y, ..., groups = NULL, subscripts) {
                 if (is.null(groups)) list(ylim= c(-1, 1) * max(abs(y[-1])))
                 else list(ylim = c(-1, 1) * max(sapply(split(y, groups[subscripts]),
                           function(x)
                           max(abs(x[-1]), na.rm = TRUE ))))
             },
             panel = panel.acfplot,
             type = if (groups) 'b' else 'h',
             aspect = "xy",
             start = 1, thin = 1,
             lag.max = NULL,
             ylab = "Autocorrelation",
             xlab = "Lag",
             main = attr(x, "title"),
             ...,
             subset = thinned.indices(x[[1]], start = start, thin = thin))
{
    if (!is.R()) {
      stop("This function is not yet available in S-PLUS")
    }
    if (groups && outer) warning("'groups=TRUE' ignored when 'outer=TRUE'")
    getAcf <- function(x, lag.max)
    {
        as.vector(acf(x, lag.max = lag.max, plot = FALSE)$acf)
    }
    if (groups || outer)
    {
        datalist <-
            lapply(x, function(x) 
                   as.data.frame(apply(as.matrix(x)[subset, ,drop=FALSE], 2,
                                       getAcf, lag.max = lag.max)))
        data <- do.call("rbind", datalist)
    }
    else
    {
        ## this is not quite valid, as we are combining multiple
        ## series, but shouldn't be too bad (FIXME: should we warn?)

        datalist <- lapply(x, function(x) as.matrix(x)[subset, ,drop=FALSE])
        data <-
            as.data.frame(apply(do.call("rbind", datalist),
                                2, getAcf, lag.max = lag.max))
    }
    form <-
        if (outer)
            as.formula(paste(paste(lapply(names(data), as.name),
                                   collapse = "+"),
                             "~ .lag | .run"))
        else
            as.formula(paste(paste(lapply(names(data), as.name),
                                   collapse = "+"),
                             "~ .lag"))
    data[[".lag"]] <- seq(length = nrow(datalist[[1]])) ## repeated
    .run <- gl(length(datalist), nrow(datalist[[1]]))
    if (groups && !outer)
        xyplot(form, data = data,
               outer = TRUE,
               groups = .run,
               prepanel = prepanel,
               panel = panel,
               type = type,
               aspect = aspect,
               xlab = xlab,
               ylab = ylab,
               main = main,
               ...)
    else
        xyplot(form, data = data,
               outer = TRUE,
               prepanel = prepanel,
               panel = panel,
               type = type,
               aspect = aspect,
               xlab = xlab,
               ylab = ylab,
               main = main,
               ...)
}


