panel.hexplom <-
    function(...)
    panel.hexbinplot(...)


hexplom <- function(x, data, ...)
{
    UseMethod("hexplom")
}





hexplom.data.frame <- 
    function (x, data = NULL, ..., groups = NULL, subset = TRUE)
{
    ocall <- sys.call(sys.parent())
    ocall[[1]] <- quote(hexplom)
    ccall <- match.call()
    if (!is.null(ccall$data))
        warning("explicit 'data' specification ignored")
    ccall$data <- list(x = x, groups = groups, subset = subset)
    ccall$x <- ~x
    ccall$groups <- groups
    ccall$subset <- subset
    ccall[[1]] <- quote(hexbin::hexplom)
    ans <- eval.parent(ccall)
    ans$call <- ocall
    ans
}

hexplom.matrix <- 
    function (x, data = NULL, ..., groups = NULL, subset = TRUE)
{
    ocall <- sys.call(sys.parent())
    ocall[[1]] <- quote(hexplom)
    ccall <- match.call()
    if (!is.null(ccall$data))
        warning("explicit 'data' specification ignored")
    ccall$data <- list(x = x, groups = groups, subset = subset)
    ccall$x <- ~x
    ccall$groups <- groups
    ccall$subset <- subset
    ccall[[1]] <- quote(hexbin::hexplom)
    ans <- eval.parent(ccall)
    ans$call <- ocall
    ans
}


hexplom.formula <-
  function(x, data = NULL, ...)
{
    ocall <- sys.call(sys.parent())
    ocall[[1]] <- quote(hexplom)
    ccall <- match.call()
    ccall[[1]] <- quote(lattice::splom)
    if (is.null(ccall$panel)) ccall$panel <- panel.hexplom
    ans <- eval.parent(ccall)
    ans$call <- ocall
    ans
}




old.hexplom.formula <-
    function(x,
             data = parent.frame(),
             auto.key = FALSE,
             aspect = 1,
             between = list(x = 0.5, y = 0.5),
             #panel = if (is.null(groups)) "panel.hexplom"
             #else "panel.superpose",
             panel = panel.hexplom,
             prepanel = NULL,
             scales = list(),
             strip = TRUE,
             groups = NULL,
             xlab = "Scatter Plot Matrix",
             xlim,
             ylab = NULL,
             ylim,
             superpanel = "panel.pairs",
             pscales = 5,
             varnames,
             drop.unused.levels = lattice.getOption("drop.unused.levels"),
             ...,
             default.scales = list(draw = FALSE, relation = "same", axs = "i"),
             subset = TRUE)
{
    ## dots <- eval(substitute(list(...)), data, parent.frame())
    dots <- list(...)

    #groups <- eval(substitute(groups), data, parent.frame())
    if(!is.null(groups))stop("groups not implemented yet")
    subset <- eval(substitute(subset), data, parent.frame())

    ## Step 1: Evaluate x, y, etc. and do some preprocessing

    ## right.name <- deparse(substitute(formula))
    ## formula <- eval(substitute(formula), data, parent.frame())
    form <-
        ## if (inherits(formula, "formula"))
        latticeParseFormula(x, data,
                            subset = subset, groups = groups,
                            multiple = FALSE,
                            outer = FALSE, subscripts = TRUE,
                            drop = drop.unused.levels)
##         else {
##             if (is.matrix(formula)) {
##                 list(left = NULL,
##                      right = as.data.frame(formula)[subset,],
##                      condition = NULL,
##                      left.name = "",
##                      right.name =  right.name,
##                      groups = groups,
##                      subscr = seq(length = nrow(formula))[subset])
##             }
##             else if (is.data.frame(formula)) {
##                 list(left = NULL,
##                      right = formula[subset,],
##                      condition = NULL,
##                      left.name = "",
##                      right.name =  right.name,
##                      groups = groups,
##                      subscr = seq(length = nrow(formula))[subset])
##             }
##             else stop("invalid formula")
##         }


    ## We need to be careful with subscripts here. It HAS to be there,
    ## and it's to be used to index x, y, z (and not only groups,
    ## unlike in xyplot etc). This means we have to subset groups as
    ## well, which is about the only use for the subscripts calculated
    ## in latticeParseFormula, after which subscripts is regenerated
    ## as a straight sequence indexing the variables

    if (!is.null(form$groups)) groups <-  form$groups[form$subscr]
    subscr <- seq(length = nrow(form$right))

    if (!is.function(panel)) panel <- eval(panel)
    if (!is.function(strip)) strip <- eval(strip)

    prepanel <-
        if (is.function(prepanel)) prepanel
        else if (is.character(prepanel)) get(prepanel)
        else eval(prepanel)

    cond <- form$condition
    number.of.cond <- length(cond)
    x <- as.data.frame(form$right)

    if (number.of.cond == 0) {
        strip <- FALSE
        cond <- list(as.factor(rep(1, nrow(x))))
        number.of.cond <- 1
    }

    if (!missing(varnames)) colnames(x) <-
        eval(substitute(varnames), data, parent.frame())

    ## create a skeleton trellis object with the
    ## less complicated components:

    #foo <- do.call(lattice:::trellis.skeleton,
    foo <- do.call(trellis.skeleton,
                   c(list(cond = cond,
                          aspect = aspect,
                          between = between,
                          panel = superpanel,
                          strip = strip,
                          xlab = xlab,
                          ylab = ylab,
                          xlab.default = "Scatter Plot Matrix"), dots))

    dots <- foo$dots # arguments not processed by trellis.skeleton
    foo <- foo$foo
    foo$call <- match.call()

    ## Step 2: Compute scales.common (leaving out limits for now)

    ## FIXME: It is not very clear exactly what effect scales is
    ## supposed to have. Not much in Trellis (probably), but there are
    ## certain components which are definitely relevant, and certain
    ## others (like log) which can be used in innovative
    ## ways. However, I'm postponing all that to later, if at all

    if (!is.list(scales)) scales <- list()

    ## some defaults for scales

#     if (is.null(scales$draw)) scales$draw <- FALSE
#     if (is.null(scales$relation)) scales$relation <- "same"
#     if (is.null(scales$axs)) scales$axs <- "i"

    scales <- updateList(default.scales, scales)
    foo <- c(foo,
             #do.call(lattice:::construct.scales, scales))
             do.call(construct.scales, scales))


    ## Step 3: Decide if limits were specified in call:

    have.xlim <- !missing(xlim)
    if (!is.null(foo$x.scales$limit)) {
        have.xlim <- TRUE
        xlim <- foo$x.scales$limit
    }
    have.ylim <- !missing(ylim)
    if (!is.null(foo$y.scales$limit)) {
        have.ylim <- TRUE
        ylim <- foo$y.scales$limit
    }
    
    ## Step 4: Decide if log scales are being used (has to be NO):

    have.xlog <- !is.logical(foo$x.scales$log) || foo$x.scales$log
    have.ylog <- !is.logical(foo$y.scales$log) || foo$y.scales$log

    ## immaterial, since scales has no effect.

#    if (have.xlog) {
#        xlog <- foo$x.scales$log
#        xbase <-
#            if (is.logical(xlog)) 10
#            else if (is.numeric(xlog)) xlog
#            else if (xlog == "e") exp(1)
#
#        x <- log(x, xbase)
#        if (have.xlim) xlim <- log(xlim, xbase)
#    }
#    if (have.ylog) {
#        ylog <- foo$y.scales$log
#        ybase <-
#            if (is.logical(ylog)) 10
#            else if (is.numeric(ylog)) ylog
#            else if (ylog == "e") exp(1)
#
#        y <- log(y, ybase)
#        if (have.ylim) ylim <- log(ylim, ybase)
#    }
    
    ## Step 5: Process cond

    cond.max.level <- unlist(lapply(cond, nlevels))

    ## id.na used only to see if any plotting is needed. Not used
    ## subsequently, unlike other functions

    id.na <- FALSE
    for (j in 1:ncol(x))
        id.na <- id.na | is.na(x[,j])
    for (var in cond)
        id.na <- id.na | is.na(var)
    if (!any(!id.na)) stop("nothing to draw")
    ## Nothing simpler ?


    ## Step 6: Evaluate layout, panel.args.common and panel.args


    foo$panel.args.common <-
        c(list(z = x,
               panel = panel,
               panel.subscripts = TRUE,
               groups = groups, # xscales = foo$x.scales, yscales =foo$y.scales,
               .aspect.ratio=aspect,
               pscales = pscales),
          dots)

    nplots <- prod(cond.max.level)
    if (nplots != prod(sapply(foo$condlevels, length))) stop("mismatch")
    foo$panel.args <- vector(mode = "list", length = nplots)


    cond.current.level <- rep(1, number.of.cond)


    for (panel.number in seq(length = nplots))
    {

        ##id <- !id.na  WHY ?
        for(i in 1:number.of.cond)
        {
            var <- cond[[i]]
            id <- if (is.shingle(var))
                ((var >=
                  levels(var)[[cond.current.level[i]]][1])
                 & (var <=
                    levels(var)[[cond.current.level[i]]][2]))
            else (as.numeric(var) == cond.current.level[i])
        }

        foo$panel.args[[panel.number]] <-
            list(subscripts = subscr[id])

        cond.current.level <-
            #lattice:::cupdate(cond.current.level, cond.max.level)
            cupdate(cond.current.level, cond.max.level)
    }


    #more.comp <- c(lattice:::limits.and.aspect(
    more.comp <- c(limits.and.aspect(
                                     lattice::prepanel.default.splom,
                                     prepanel = prepanel, 
                                     have.xlim = have.xlim, xlim = xlim, 
                                     have.ylim = have.ylim, ylim = ylim, 
                                     x.relation = foo$x.scales$relation,
                                     y.relation = foo$y.scales$relation,
                                     panel.args.common = foo$panel.args.common,
                                     panel.args = foo$panel.args,
                                     aspect = aspect,
                                     nplots = nplots,
                                     x.axs = foo$x.scales$axs,
                                     y.axs = foo$y.scales$axs),
                  #lattice::: cond.orders(foo))
                  cond.orders(foo))
    foo[names(more.comp)] <- more.comp



    if (is.null(foo$legend) && !is.null(groups) &&
        (is.list(auto.key) || (is.logical(auto.key) && auto.key)))
    {
        foo$legend <-
            list(list(fun = "drawSimpleKey",
                      args =
                      updateList(list(text = levels(as.factor(groups)),
                                      points = TRUE,
                                      rectangles = FALSE,
                                      lines = FALSE), 
                                 if (is.list(auto.key)) auto.key else list())))
        foo$legend[[1]]$x <- foo$legend[[1]]$args$x
        foo$legend[[1]]$y <- foo$legend[[1]]$args$y
        foo$legend[[1]]$corner <- foo$legend[[1]]$args$corner

        names(foo$legend) <- 
            if (any(c("x", "y", "corner") %in% names(foo$legend[[1]]$args)))
                "inside"
            else
                "top"
        if (!is.null(foo$legend[[1]]$args$space))
            names(foo$legend) <- foo$legend[[1]]$args$space
    }

    class(foo) <- "trellis"
    foo
}
