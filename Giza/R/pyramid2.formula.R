pyramid2.formula <-
function (x, data = NULL, allow.multiple = is.null(groups) || 
    outer, outer = FALSE, auto.key = FALSE, aspect = "fill", 
    panel = panel, prepanel = NULL, 
    scales = list(), strip = TRUE, groups = NULL, xlab, xlim, 
    ylab, ylim, box.ratio = 1, horizontal = NULL, drop.unused.levels = lattice.getOption("drop.unused.levels"), 
    ..., lattice.options = NULL, default.scales = if (horizontal) list(y = list(tck = 0, 
        alternating = FALSE, rot = 0)) else list(x = list(tck = 0, 
        alternating = FALSE)), default.prepanel = prepanel, 
    subscripts = !is.null(groups), subset = TRUE) 
{

    formula <- x
    dots <- list(...)
    groups <- eval(substitute(groups), data, environment(formula))
    subset <- eval(substitute(subset), data, environment(formula))
    if (!is.null(lattice.options)) { 
        oopt <- lattice.options(lattice.options)
        on.exit(lattice.options(oopt), add = TRUE)
    }
    is.standard.barchart <- is.character(panel) && panel == "panel.barchart2"
    form <- lattice:::latticeParseFormula(formula, data, subset = subset, 
        groups = groups, multiple = allow.multiple, outer = outer, 
        subscripts = TRUE, drop = drop.unused.levels)
    groups <- form$groups
    if (!is.function(panel)) panel <- eval(panel)
    if (!is.function(strip)) strip <- eval(strip)
    if ("subscripts" %in% names(formals(panel))) subscripts <- TRUE
    if (subscripts) subscr <- form$subscr
    cond <- form$condition
    x <- form$right
    y <- form$left
    if (is.null(y)) {
        y <- rep(if (is.null(names(x))) "" else names(x), length.out = length(x))
        y <- factor(y, levels = unique(y))
    }
    if (length(cond) == 0) {
        strip <- FALSE
        cond <- list(gl(1, length(x)))
    }
    if (is.null(horizontal)) {
        horizontal <- if ((is.factor(x) || is.shingle(x) || is.character(x)) && 
            is.numeric(y)) 
            FALSE
        else TRUE
    }
    if (horizontal) {
        y <- as.factorOrShingle(y)
        is.f.y <- is.factor(y)
        num.l.y <- nlevels(y)
        if (missing(xlab)) 
            xlab <- form$right.name
        if (missing(ylab)) 
            ylab <- if (is.f.y) 
                NULL
            else form$left.name
    } else {
        x <- as.factorOrShingle(x)
        is.f.x <- is.factor(x)
        num.l.x <- nlevels(x)
        if (missing(ylab)) 
            ylab <- form$left.name
        if (missing(xlab)) 
            xlab <- if (is.f.x) 
                NULL
            else form$right.name
    }
    trellis.skeleton <- lattice:::trellis.skeleton 
    foo <- do.call("trellis.skeleton", c(list(formula = formula, 
        cond = cond, aspect = aspect, strip = strip, panel = panel, 
        xlab = xlab, ylab = ylab, xlab.default = form$right.name, 
        ylab.default = form$left.name, lattice.options = lattice.options), 
        dots))
  
    dots <- foo$dots
    foo <- foo$foo
    foo$call <- sys.call(sys.parent())
    foo$call[[1]] <- quote(bwplot)  
    if (is.character(scales)) scales <- list(relation = scales)
    scales <- lattice:::updateList(default.scales, scales)
    construct.scales <- lattice:::construct.scales
    foo <- c(foo, do.call("construct.scales", scales))
    have.xlim <- !missing(xlim)
    if (!is.null(foo$x.scales$limits)) {
        have.xlim <- TRUE
        xlim <- foo$x.scales$limits
    }
    have.ylim <- !missing(ylim)
    if (!is.null(foo$y.scales$limits)) {
        have.ylim <- TRUE
        ylim <- foo$y.scales$limits
    }
    have.xlog <- !is.logical(foo$x.scales$log) || foo$x.scales$log
    have.ylog <- !is.logical(foo$y.scales$log) || foo$y.scales$log
    if (have.xlog) {
        xlog <- foo$x.scales$log
        xbase <- if (is.logical(xlog)) 
            10
        else if (is.numeric(xlog)) 
            xlog
        else if (xlog == "e") 
            exp(1)
        x <- log(x, xbase)
        if (have.xlim) 
            xlim <- lattice:::logLimits(xlim, xbase)
    }
    if (have.ylog) {
        ylog <- foo$y.scales$log
        ybase <- if (is.logical(ylog)) 
            10
        else if (is.numeric(ylog)) 
            ylog
        else if (ylog == "e") 
            exp(1)
        y <- log(y, ybase)
        if (have.ylim) 
            ylim <- lattice:::logLimits(ylim, ybase)
    }
    cond.max.level <- unlist(lapply(cond, nlevels))
    foo$panel.args.common <- dots
    foo$panel.args.common$box.ratio <- box.ratio
    foo$panel.args.common$horizontal <- horizontal
    if (subscripts) foo$panel.args.common$groups <- groups
    if (horizontal) {
        if (!is.f.y) 
            foo$panel.args.common$nlevels <- num.l.y
    } else {
        if (!is.f.x) foo$panel.args.common$nlevels <- num.l.x
    }
    npackets <- prod(cond.max.level)
    if (npackets != prod(sapply(foo$condlevels, length))) stop("mismatch in number of packets")
    foo$panel.args <- vector(mode = "list", length = npackets)
    foo$packet.sizes <- numeric(npackets)
    if (npackets > 1) {
        dim(foo$packet.sizes) <- sapply(foo$condlevels, length)
        dimnames(foo$packet.sizes) <- lapply(foo$condlevels, 
            as.character)
    }
 
  cond.current.level <- rep(1, length(cond))
    for (packet.number in seq_len(npackets)) {
        id <- lattice:::compute.packet(cond, cond.current.level)
        foo$packet.sizes[packet.number] <- sum(id)
        if (horizontal) {
            if (is.f.y) {print("####### you have to fix something here first if you want to run it this way #######")
                foo$panel.args[[packet.number]] <- list(x = x[id], 
                  y = y[id])
                if (subscripts) 
                  foo$panel.args[[packet.number]]$subscripts <- subscr[id]
            } else {
                panel.x <- numeric(0)
                panel.y <- numeric(0)
                if (subscripts) 
                  panel.subscr <- numeric(0)
                for (k in seq_len(num.l.y)) {
                  tid <- id & (y >= levels(y)[[k]][1]) & (y <= 
                    levels(y)[[k]][2])
                  panel.x <- c(panel.x, x[tid])
                  panel.y <- c(panel.y, rep(k, length(tid[tid])))
                  if (subscripts) 
                    panel.subscr <- c(panel.subscr, subscr[tid])
                }
                foo$panel.args[[packet.number]] <- list(x = panel.x, 
                  y = panel.y)
                if (subscripts) 
                  foo$panel.args[[packet.number]]$subscripts <- panel.subscr
            }
        } else {
            if (is.f.x) {
                foo$panel.args[[packet.number]] <- list(x = x[id], 
                  y = y[id])
                if (subscripts) 
                  foo$panel.args[[packet.number]]$subscripts <- subscr[id]
            }
            else {
                panel.x <- numeric(0)
                panel.y <- numeric(0)
                if (subscripts) 
                  panel.subscr <- numeric(0)
                for (k in seq_len(num.l.x)) {
                  tid <- id & (x >= levels(x)[[k]][1]) & (x <= 
                    levels(x)[[k]][2])
                  panel.y <- c(panel.y, y[tid])
                  panel.x <- c(panel.x, rep(k, length(tid[tid])))
                  if (subscripts) 
                    panel.subscr <- c(panel.subscr, subscr[tid])
                }
                foo$panel.args[[packet.number]] <- list(x = panel.x, 
                  y = panel.y)
                if (subscripts) 
                  foo$panel.args[[packet.number]]$subscripts <- panel.subscr
            }
        }
        cond.current.level <- lattice:::cupdate(cond.current.level, cond.max.level)
    }
    more.comp <- c(lattice:::limits.and.aspect(prepanel.default = prepanel.default.bwplot2, prepanel = prepanel, 
        have.xlim = have.xlim, xlim = xlim, have.ylim = have.ylim, 
        ylim = ylim, x.relation = foo$x.scales$relation, y.relation = foo$y.scales$relation, 
        panel.args.common = foo$panel.args.common, panel.args = foo$panel.args, 
        aspect = aspect, npackets = npackets, x.axs = foo$x.scales$axs, 
        y.axs = foo$y.scales$axs), lattice:::cond.orders(foo))
    foo[names(more.comp)] <- more.comp
    if (is.null(foo$legend) && lattice:::needAutoKey(auto.key, groups)) {
        foo$legend <- list(list(fun = "drawSimpleKey", args = lattice:::updateList(list(text = levels(as.factor(groups)), 
            points = if (is.standard.barchart) FALSE else TRUE, 
            rectangles = if (is.standard.barchart) TRUE else FALSE, 
            lines = FALSE), if (is.list(auto.key)) auto.key else list())))
        foo$legend[[1]]$x <- foo$legend[[1]]$args$x
        foo$legend[[1]]$y <- foo$legend[[1]]$args$y
        foo$legend[[1]]$corner <- foo$legend[[1]]$args$corner
        names(foo$legend) <- if (any(c("x", "y", "corner") %in% 
            names(foo$legend[[1]]$args))) 
            "inside"
        else "top"
        if (!is.null(foo$legend[[1]]$args$space)) 
            names(foo$legend) <- foo$legend[[1]]$args$space
    }
    class(foo) <- "trellis"
    foo
}
