## copies of unexported functions from R packages
## needed for R-devel leading to R-3.1.0

## used in t.trellis.R
## > lattice:::t.trellis
lattice.t.trellis <-
function (x)
{
    stopifnot(length(dim(x)) == 2)
    update(x, perm.cond = rev(x$perm.cond))
}
## <bytecode: 0x076f8c8c>
## <environment: namespace:lattice>



## used in grid.axis.hh.R
## > grid:::make.yaxis.major
grid.make.yaxis.major <-
function (at, main)
{
    if (main)
        x <- c(0, 0)
    else x <- c(1, 1)
    linesGrob(unit(x, "npc"), unit(c(min(at), max(at)), "native"),
        name = "major")
}
## <bytecode: 0x076f6d68>
## <environment: namespace:grid>


## > grid:::make.yaxis.ticks
grid.make.yaxis.ticks <-
function (at, main)
{
    if (main) {
        tick.x0 <- unit(0, "npc")
        tick.x1 <- unit(-0.5, "lines")
    }
    else {
        tick.x0 <- unit(1, "npc")
        tick.x1 <- unit(1, "npc") + unit(0.5, "lines")
    }
    segmentsGrob(tick.x0, unit(at, "native"), tick.x1, unit(at,
        "native"), name = "ticks")
}
## <bytecode: 0x076f4098>
## <environment: namespace:grid>


## > grid:::make.xaxis.major
grid.make.xaxis.major <-
function (at, main)
{
    if (main)
        y <- c(0, 0)
    else y <- c(1, 1)
    linesGrob(unit(c(min(at), max(at)), "native"), unit(y, "npc"),
        name = "major")
}
## <bytecode: 0x076ef630>
## <environment: namespace:grid>


## > grid:::make.xaxis.ticks
grid.make.xaxis.ticks <-
function (at, main)
{
    if (main) {
        tick.y0 <- unit(0, "npc")
        tick.y1 <- unit(-0.5, "lines")
    }
    else {
        tick.y0 <- unit(1, "npc")
        tick.y1 <- unit(1, "npc") + unit(0.5, "lines")
    }
    segmentsGrob(unit(at, "native"), tick.y0, unit(at, "native"),
        tick.y1, name = "ticks")
}
## <bytecode: 0x076ed77c>
## <environment: namespace:grid>



## used in glht.mmc.R and mcalinfct.R
## > multcomp:::meanslinfct
multcomp.meanslinfct <-
function (model, focus, mmm.data = model$model, formula.in = terms(model),
    contrasts.arg = NULL)
{
    mmm.factor <- sapply(mmm.data, inherits, "factor")
    mmm.levels <- lapply(mmm.data[mmm.factor], levels)
    mmm.rows <- sapply(mmm.levels, length)
    n.mmm.rows <- prod(mmm.rows)
    mmm.new <- mmm.data[1:n.mmm.rows, ]
    mmm.factor.names <- names(mmm.data)[mmm.factor]
    mmm.rows.forward <- cumprod(mmm.rows)
    mmm.rows.forward.prev <- c(1, mmm.rows.forward)
    names(mmm.rows.forward.prev) <- c(names(mmm.rows.forward),
        "all")
    for (i in mmm.factor.names) mmm.new[[i]] <- gl(mmm.rows[i],
        mmm.rows.forward.prev[i], n.mmm.rows, labels = mmm.levels[[i]])
    mmm.numeric.names <- names(mmm.data)[!mmm.factor]
    for (i in mmm.numeric.names) mmm.new[[i]][] <- mean(mmm.data[[i]])
    none.data <- model.matrix(formula.in, data = mmm.new, contrasts.arg = contrasts.arg)
    none.linfct <- aggregate(none.data, by = mmm.new[focus],
        FUN = mean)[, -1]
    rownames(none.linfct) <- levels(mmm.new[[focus]])
    data.matrix(none.linfct)
}
## <environment: namespace:multcomp>


## used in ae.dotplot7.R
## > lattice:::lattice.setStatus
lattice.lattice.setStatus <-
function (..., prefix = NULL, clean.first = FALSE)
{
    dots <- list(...)
    if (is.null(names(dots)) && length(dots) == 1 && is.list(dots[[1]]))
        dots <- dots[[1]]
    if (length(dots) == 0)
        return()
    .LatticeEnv <- get(".LatticeEnv", envir=environment(lattice.options)) ## danger here
    lattice.status <- if (clean.first)
        list()
    else get("lattice.status", envir = .LatticeEnv)
    if (is.null(prefix))
        lattice.status[names(dots)] <- dots
    else lattice.status[[prefix]][names(dots)] <- dots
    assign("lattice.status", lattice.status, envir = .LatticeEnv)
    invisible()
}
## <bytecode: 0x103970ae0>
## <environment: namespace:lattice>


## > lattice:::hist.constructor
lattice.hist.constructor <-
function (x, breaks, include.lowest = TRUE, right = TRUE, ...)
{
    if (is.numeric(breaks) && length(breaks) > 1)
        hist(as.numeric(x), breaks = breaks, plot = FALSE, include.lowest = include.lowest,
            right = right)
    else hist(as.numeric(x), breaks = breaks, right = right,
        plot = FALSE)
}
## <bytecode: 0x105f79378>
## <environment: namespace:lattice>


## used in panel.axis.right.R
## lattice:::extend.limits
## lattice:::chooseFace
## lattice:::lattice.getStatus

## > lattice:::extend.limits
lattice.extend.limits <-
function (lim, length = 1, axs = "r", prop = if (axs == "i") 0 else lattice.getOption("axis.padding")$numeric)
{
    if (all(is.na(lim)))
        NA_real_
    else if (is.character(lim)) {
        c(1, length(lim)) + c(-1, 1) * if (axs == "i")
            0.5
        else lattice.getOption("axis.padding")$factor
    }
    else if (length(lim) == 2) {
        if (lim[1] > lim[2]) {
            ccall <- match.call()
            ccall$lim <- rev(lim)
            ans <- eval.parent(ccall)
            return(rev(ans))
        }
        if (!missing(length) && !missing(prop))
            stop("'length' and 'prop' cannot both be specified")
        if (length <= 0)
            stop("'length' must be positive")
        if (!missing(length)) {
            prop <- (as.numeric(length) - as.numeric(diff(lim)))/(2 *
                as.numeric(diff(lim)))
        }
        if (lim[1] == lim[2])
            lim + 0.5 * c(-length, length)
        else {
            d <- diff(as.numeric(lim))
            lim + prop * d * c(-1, 1)
        }
    }
    else {
        print(lim)
        stop("improper length of 'lim'")
    }
}
## <bytecode: 0x1060db350>
## <environment: namespace:lattice>


## > lattice:::chooseFace
lattice.chooseFace <-
function (fontface = NULL, font = 1)
{
    if (is.null(fontface))
        font
    else fontface
}
## <bytecode: 0x102add8e8>
## <environment: namespace:lattice>


## > lattice:::lattice.getStatus
lattice.lattice.getStatus <-
function (name, prefix = NULL)
{
    .LatticeEnv <- get(".LatticeEnv", envir=environment(lattice.options)) ## danger here
    if (is.null(prefix))
        get("lattice.status", envir = .LatticeEnv)[[name]]
    else get("lattice.status", envir = .LatticeEnv)[[prefix]][[name]]
}
## <bytecode: 0x10620ff58>
## <environment: namespace:lattice>

