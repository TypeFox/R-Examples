# the functions in this file are verbatim copies from those in package
# lattice, http://cran.r-project.org/src/contrib/lattice_0.20-29.tar.gz
# copied on Aug 8, 2014, by Edzer Pebesma.

# reason for copying is that hexbin 1.26-3 generates
# the following NOTE on CRAN:
# 
# checking dependencies in R code ... NOTE
# Unexported objects imported by ':::' calls:
# lattice:::cond.orders lattice:::construct.scales
# lattice:::cupdate lattice:::limits.and.aspect
# lattice:::trellis.skeleton
# See the note in ?::: about the use of this operator.
# See the information on DESCRIPTION files in the chapter Creating R
# packages of the Writing R Extensions manual.

# the files in lattice carry the following copyright notice:

### Copyright (C) 2001-2006  Deepayan Sarkar <Deepayan.Sarkar@R-project.org>
### Copyright (C) 2001-2005  Saikat DebRoy <saikat@stat.wisc.edu>
###
### This file is part of the lattice package for R.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
### MA 02110-1301, USA

construct.legend <-
    function(legend = NULL, key = NULL, fun = "draw.key")
{
    if (is.null(legend) && is.null(key)) return(NULL)
    if (is.null(legend)) legend <- list()
    if (!is.null(key))
    {
        space <- key$space
        x <- y <- corner <- NULL

        if (is.null(space))
        {
            if (any(c("x", "y", "corner") %in% names(key)))
            {
                stopifnot(is.null(x) || (length(x) == 1 && x >= 0 && x <= 1))
                stopifnot(is.null(y) || (length(y) == 1 && y >= 0 && y <= 1))
                stopifnot(is.null(corner) ||
                          (length(corner) == 2 &&
                           all(corner %in% c(0, 1))))
                space <- "inside"
                x <- key$x
                y <- key$y
                corner <- key$corner
                ## check for valid values
            }
            else
                space <- "top"
        }
        if (space != "inside" && space %in% names(legend))
            stop(gettextf("component '%s' duplicated in key and legend", space))

        key.legend <- list(fun = fun, args = list(key = key, draw = FALSE))
        key.legend$x <- x
        key.legend$y <- y
        key.legend$corner <- corner

        legend <- c(list(key.legend), legend)
        names(legend)[1] <- space
    }
    legend
}

extend.limits <-
    function(lim, length = 1, axs = "r",
             prop =
             if (axs == "i") 0
             else lattice.getOption("axis.padding")$numeric)
{
    ## if (!is.numeric(lim)) NA
    if (all(is.na(lim))) NA_real_ # or lim?
    else if (is.character(lim) )
    {
        c(1, length(lim)) + c(-1, 1) * if (axs == "i") 0.5 else lattice.getOption("axis.padding")$factor
    }
    else if (length(lim) == 2)
    {
        if (lim[1] > lim[2])
        {
            ccall <- match.call()
            ccall$lim <- rev(lim)
            ans <- eval.parent(ccall)
            return (rev(ans))
        }
        if (!missing(length) && !missing(prop))
            stop("'length' and 'prop' cannot both be specified")
        if (length <= 0) stop("'length' must be positive")
        if (!missing(length))
        {
            prop <- (as.numeric(length) - as.numeric(diff(lim))) / (2 * as.numeric(diff(lim)))
        }
        if (lim[1]==lim[2]) lim + 0.5 * c(-length,length)
        else
        {
            d <- diff(as.numeric(lim))
            lim + prop * d * c(-1,1)
        }
    }
    else
    {
        print(lim)
        stop("improper length of 'lim'")
    }
}

limitsFromLimitlist <-
    function(have.lim,
             lim,
             relation,
             limitlist,
             used.at,
             numlimitlist,
             axs,
             npackets)
    ## have.lim: logical, whether xlim/ylim was explicitly specified
    ## lim: the specified limit if have.lim = TRUE
    ## relation: same/free/sliced
    ## limitlist: list of limits from prepanel calculations, one for each panel
    ## numlimitlist: (optional) numeric locations for factors (lim
    ##                will be levels including unused ones)
    ## axs: "r", "i" etc, passed on to extend.limits

    ## return value depends on relation. (See limits.and.aspect below,
    ## where this is used, for partial enlightenment.)
{

    if (relation == "same")
    {
        ## The problem here is that we need to figure out the overall
        ## limit required from the limits of each panel. This could be
        ## a problem for two reasons. First, some panels could have no
        ## data in them, in which case the corresponding limits would
        ## be NA. Secondly, the limits could be either numeric or
        ## character vectors (the latter for factors). When relation =
        ## same, the type should be same across panels. When numeric,
        ## we just take range, leaving out NAs. But what about
        ## factors?  Is it OK to assume that all the non-NA vectors
        ## would be exactly the same ? They should be, since levels(x)
        ## would not change even if not all levels are
        ## represented. So, I'm just taking unique of all the vectors
        ## concatenated, excluding NA's

        ## Additional complication: Need to preserve class of limits,
        ## to be used later in tick location/label calculation. Not a
        ## problem in other cases, but here unlist-ing loses the
        ## class.


        #if (!have.lim)
        ## always calculate the limits from prepanel first:

        ## should check that all classes are the same. How ? What
        ## about NA's ? Arrgh!

        ## to handle NA's, how about:

        all.na <- unlist(lapply(limitlist, function(x) all(is.na(x))))
        class.lim <- ## retain non-NA limitlists only
            lapply(limitlist[!all.na], class)
        ## class.lim is a list now, may be length 0
        limits <- unlist(limitlist) ## loses the class attribute

        ## if (length(limits) > 0)
        if (sum(!is.na(limits)) > 0)
        {
            if (is.character(limits))
            {
                limits <- unique(limits[!is.na(limits)])
                slicelen <- diff(extend.limits(limits, axs = axs))
            }
            else ## if (is.numeric(limits)) # or dates etc
            {
                limits <-
                    extend.limits(range(as.numeric(limits), finite = TRUE),
                                  axs = axs)
                slicelen <- diff(range(limits, finite = TRUE))
            }

            ## hopefully put back appropriate class of limits:
            ## FIXME: date changes may have messed this up

            if (length(class.lim) > 0)
                class(limits) <-
                    if (all(class.lim[[1]] == "integer"))
                        "numeric" else class.lim[[1]]

            ## (have to handle "integer" specially, since variable
            ## specifications like 1:10 are rather common, and
            ## class() <- "integer" would turn the limits into
            ## integers)
        }
        else
        {
            limits <- c(0,1)
            slicelen <- 1
        }

        if (have.lim)
        {
            if (is.list(lim))
                stop("limits cannot be a list when relation = same")
            old.limits <- limits
            limits <- lim
            ## lim overrides prepanel except NAs
            if (!is.character(limits) && !is.character(old.limits)) {
                limits[is.na(limits)] <- old.limits[is.na(limits)]
            }
            slicelen <-
                ## this no longer works for dates (R 2.6)
##                 if (is.numeric(lim)) diff(range(lim))
##                 else length(lim) + 2
                if (is.character(limits)) length(limits) + 2
                else diff(range(as.numeric(limits)))
        }
        ans <- list(limits = limits, slicelen = slicelen)
    }
    else if (relation == "sliced")
    {
        if (have.lim)
        {
            if (is.list(lim))
            {
                limits <- rep(lim, length.out = npackets)
            }
            else warning("Explicitly specified limits ignored")
        }
        slicelen <- limitlist
        for (i in seq_along(limitlist))
        {
            slicelen[[i]] <-
                ## if (is.numeric(limitlist[[i]]))
                if (!is.character(limitlist[[i]]))
                {
                    if (any(is.finite(limitlist[[i]])))
                        ## range unnecessary, but...
                        diff(range(as.numeric(limitlist[[i]]), finite = TRUE))
                    else NA_real_
                }
                else if (!any(is.na(numlimitlist[[i]])))
                    diff(range(as.numeric(numlimitlist[[i]])))
                else NA_real_
        }
        slicelen <-
            (if (axs == "i") 1 else 1 + 2 * lattice.getOption("axis.padding")$numeric) *
                max(unlist(slicelen), na.rm = TRUE)
        for (i in seq_along(limitlist))
        {
            if (is.numeric(limitlist[[i]]))
                limitlist[[i]] <-
                    extend.limits(limitlist[[i]], length = slicelen)
        }
        for (i in seq_along(numlimitlist))
        {
            if (!all(is.na(numlimitlist[[i]])))
                numlimitlist[[i]] <-
                    extend.limits(as.numeric(numlimitlist[[i]]), length = slicelen)
        }
        ans <-
            list(limits = limitlist,
                 used.at = used.at,
                 numlimitlist = numlimitlist,
                 slicelen = slicelen)
    }
    else if (relation == "free")
    {
        if (have.lim)
        {
            ## This is the only situation where limits can be a list
            ## (doesn't make sense when relation="same", ignored when
            ## relation="sliced").  Even if limits is not a list (but
            ## is specified), it will be treated as a list, and
            ## repeated as necessary (see further comments below).

            if (!is.list(lim)) lim <- list(lim)

            ## There's a subtle consideration here.  It is possible
            ## for some panels to have nothing in them (or only NA's).
            ## Such panels usually have their prepanel functions
            ## return NA.  When 'limits' is specified as a list, this
            ## will be interpreted as the limit specification for the
            ## non-empty panels only (this is an arbitrary choice, but
            ## it usually makes more sense, even though it's less
            ## general than the other choice).

            ## which ones are non-NA?
            id <- which(sapply(limitlist, function(x) !all(is.na(x))))

            ## replace these with the limits supplied, except if the
            ## supplied limits are NULL, in which case retain limits
            ## calculated by prepanel.

            old.limitlist <- limitlist
            limitlist[id] <- lim
            which.null <- sapply(limitlist, is.null)
            limitlist[which.null] <- old.limitlist[which.null]

            ## lim overrides prepanel except NAs
            for (i in seq_along(limitlist))
            {
                if (!is.character(limitlist[[i]]) &&
                    !is.character(old.limitlist[[i]]))
                {
                    isna <- is.na(limitlist[[i]])
                    limitlist[[i]][isna] <- old.limitlist[[i]][isna]
                }
            }
        }
        for (i in seq_along(limitlist))
        {
            if (!all(is.na(limitlist[[i]])) && !is.character(limitlist[[i]])) 
                limitlist[[i]] <- ## preserves class
                    extend.limits(limitlist[[i]], axs = axs)
            ## o.w., keep it as it is
        }
        slicelen <- numeric(length(limitlist))
        for (i in seq_along(limitlist))
            slicelen[i] <-
                if (!is.character(limitlist[[i]]))
                    diff(range(as.numeric(limitlist[[i]])))
                else if (!any(is.na(numlimitlist[[i]])))
                    diff(range(numlimitlist[[i]]))
                else NA_real_
        ans <-
            list(limits = limitlist,
                 used.at = used.at,
                 numlimitlist = numlimitlist,
                 slicelen = slicelen)
    }
    ans
}

complete_names <- function(x, template, allow.invalid = FALSE)
{
    pid <- pmatch(names(x), names(template), duplicates.ok = TRUE)
    if (allow.invalid) {
        x <- x[!is.na(pid)]
        pid <- pid[!is.na(pid)]
    } else {
        if (any(is.na(pid)))
            warning("Invalid or ambiguous component names: ",
                     paste(names(x)[which(is.na(pid))], collapse = ", ") )
    }
    if (any(duplicated(pid))) stop("Multiple matches to component name")
    names(x) <- names(template)[pid]
    x
}

getFunctionOrName <- function(FUN)
     ## Try lattice namespace first? Does that happen automatically?
{
    if (is.function(FUN)) FUN
    else if (is.character(FUN)) get(FUN)
    else eval(FUN)
}

trellis.skeleton <-
    function(formula = NULL,
             cond,
             aspect = default.args$aspect, # argument in xyplot
             as.table = default.args$as.table,
             between = default.args$between,
             key = NULL,
             legend = NULL,
             page = default.args$page,
             main = default.args$main,
             sub = default.args$sub,
             par.strip.text = default.args$par.strip.text,
             layout = default.args$layout,
             skip = default.args$skip,
             strip = default.args$strip.default, # argument in xyplot
             strip.left = FALSE,
             xlab.default = NULL,
             ylab.default = NULL,
             xlab = NULL, # argument in xyplot
             ylab = NULL, # argument in xyplot
             xlab.top = NULL,
             ylab.right = NULL,

             panel,       # argument in xyplot

             xscale.components = default.args$xscale.components,
             yscale.components = default.args$yscale.components,
             axis = default.args$axis,

             subscripts = TRUE, # ignored, for reasons given above

             index.cond = NULL,
             perm.cond = NULL,
             ...,
             par.settings = NULL,
             plot.args = NULL,
             lattice.options = NULL)
{
    default.args <- lattice.getOption("default.args")
    if (is.null(skip)) skip <- FALSE
    foo <-
        list(formula = formula,
             as.table = as.table,
             aspect.fill = (aspect == "fill"),
             ## key = key,
             legend = construct.legend(legend = legend, key = key),
             panel = panel,
             page = page,
             layout = layout,
             skip = skip,
             strip = if (is.logical(strip) && strip) "strip.default"
             else strip,
             strip.left = if (is.logical(strip.left) && strip.left) strip.custom(horizontal = FALSE)
             else strip.left,
             xscale.components = xscale.components,
             yscale.components = yscale.components,
             axis = axis,
             xlab = xlab,
             ylab = ylab,
             xlab.default = xlab.default,
             ylab.default = ylab.default,
             xlab.top = xlab.top,
             ylab.right = ylab.right,
             main = main,
             sub = sub,
             x.between = 0,
             y.between = 0,
             par.settings = par.settings,
             plot.args = plot.args,
             lattice.options = lattice.options,
             par.strip.text = par.strip.text,
             index.cond = index.cond,
             perm.cond = perm.cond)

    if (!is.null(between$x)) foo$x.between <- between$x
    if (!is.null(between$y)) foo$y.between <- between$y

    foo$condlevels <- lapply(cond, levels)

    list(foo = foo, dots = list(...))
}








cond.orders <- function(foo, ...)
    ## function to determine order of panels within a cond. variable
    ## foo: trellis object-to-be

    ## calculate actual values for index.cond and perm.cond.
    ## index.cond can be a function, in which case it would be used to
    ## determing order of levels within conditioning variables

    ## Question: should these be determined at run-time? Wouldn't be
    ## impossible, but has the disadvantage that looking at the
    ## trellis object will be totally uninformative in the default
    ## case (when both would be NULL). In a sense, this is fine, since
    ## having index.cond be a function is similar to having a prepanel
    ## function. After all, the results depend only on the panel
    ## contents, and those cannot be changed via update.

{

    ## the following to be used for changing order of conditioning
    ## variables and indexing their levels. The object foo already has
    ## components index.cond and perm.cond as whatever was passed to
    ## the original function call. If these are NULL, suitable
    ## defaults need to be computed. If foo$index.cond is a function,
    ## index.cond has to be computed appropriately.

    index.cond <-
        vector(mode = "list",
               length = length(foo$condlevels))

    for (i in seq_along(foo$condlevels))
        index.cond[[i]] <- seq_along(foo$condlevels[[i]])
    perm.cond <- seq_len(length(foo$condlevels))

    if (!is.null(foo$perm.cond))
    {
        if (all(sort(foo$perm.cond) == perm.cond))
            perm.cond <- foo$perm.cond
        else  stop("Invalid value of perm.cond")
    }
    if (!is.null(foo$index.cond))
    {
        if (is.list(foo$index.cond) && length(foo$index.cond) == length(index.cond))
        {
            for (i in seq_along(foo$condlevels))
                index.cond[[i]] <- index.cond[[i]][foo$index.cond[[i]]]
        }
        else if (is.function(foo$index.cond))
        {
            FUN <- foo$index.cond
            nplots <- length(foo$panel.args)
            panel.order <- numeric(nplots)
            for (count in seq_len(nplots))
            {
                if (is.list(foo$panel.args[[count]]))
                {
                    pargs <- c(foo$panel.args.common, foo$panel.args[[count]], list(...))
                    prenames <- names(formals(FUN))
                    if (!("..." %in% prenames)) pargs <- pargs[intersect(names(pargs), prenames)]
                    panel.order[count] <- do.call("FUN", pargs)
                }
                else  ## this happens for empty panels
                {
                    is.na(panel.order) <- count # panel.order[count] <- NA
                }
            }
            dim(panel.order) <- sapply(foo$condlevels, length)
            for (i in seq_along(foo$condlevels))
                index.cond[[i]] <-
                    order(apply(panel.order, i, mean, na.rm = TRUE))
        }
        else stop("Invalid value of index.cond")
    }
    list(index.cond = index.cond, perm.cond = perm.cond)
}


construct.scales <-
    function(draw = TRUE, axs = "r", tck = 1, tick.number = 5,
             at = FALSE, labels = FALSE, log = FALSE,
             alternating = TRUE, relation = "same",
             abbreviate = FALSE, minlength = 4,
             limits = NULL, format = NULL,
             equispaced.log = TRUE,

             lty = FALSE, lwd = FALSE, cex = FALSE, rot = FALSE,
             col = FALSE, col.line = col, alpha = FALSE, alpha.line = alpha,
             font = FALSE, fontfamily = FALSE, fontface = FALSE, lineheight = FALSE,

             ...,  ## NOTE: ... is currently ignored
             x = NULL, y = NULL)
{
    ## top-level values
    x.scales <- y.scales <-
        list(draw = draw, axs = axs, tck = tck, tick.number = tick.number,
             at = at, labels = labels, log = log,
             alternating = alternating, relation = relation,
             abbreviate = abbreviate, minlength = minlength,
             limits = limits, format = format, equispaced.log = equispaced.log,
             lty = lty, lwd = lwd, cex = cex, rot = rot,
             col = col, col.line = col.line, alpha = alpha, alpha.line = alpha.line,
             font = font, fontfamily = fontfamily, fontface = fontface, lineheight = lineheight)
    ## override by component-specific values
    if (!is.null(x))
    {
        if (is.character(x)) x <- list(relation = x)
        x <- complete_names(x, x.scales)
        x.scales[names(x)] <- x
    }
    if (!is.null(y))
    {
        if (is.character(y)) y <- list(relation = y)
        y <- complete_names(y, y.scales)
        y.scales[names(y)] <- y
    }
    if (is.logical(x.scales$alternating))
        x.scales$alternating <-
            if (x.scales$alternating) c(1,2)
            else 1
    if (is.logical(y.scales$alternating))
        y.scales$alternating <-
            if (y.scales$alternating) c(1,2)
            else 1
    for (nm in c("tck", "cex", "rot")) {
        x.scales[[nm]] <- rep(x.scales[[nm]], length.out = 2)
        y.scales[[nm]] <- rep(y.scales[[nm]], length.out = 2)
    }
    if (x.scales$relation == "same" && (is.list(x.scales$at) || is.list(x.scales$labels)))
        stop("the 'at' and 'labels' components of 'scales' may not be lists when 'relation = \"same\"'")
    if (y.scales$relation == "same" && (is.list(y.scales$at) || is.list(y.scales$labels)))
        stop("the 'at' and 'labels' components of 'scales' may not be lists when 'relation = \"same\"'")
    list(x.scales = x.scales, y.scales = y.scales)
}

cupdate <- function(index, maxim)
{

    ## This unexported function is used to handle arbitrary number of
    ## conditioning variables : every time it is called, it increments
    ## the "current" level of the conditioning variables suitably,
    ## i.e., it tries to increment the level of the 1st conditining
    ## variable (the one which varies fastest along panel order) and
    ## if it happens to be at its maximum (last) value, it sets it to
    ## the first value AND increments the "current" level of the 2nd
    ## (next) conditioning variable recursively.

    if(length(index)!=length(maxim)||length(maxim)<=0)
        stop("Inappropriate arguments")
    index[1] <- index[1] + 1
    if (index[1] > maxim[1] && length(maxim) > 1)
        c(1, cupdate(index[-1], maxim[-1]))
    else index
}
limits.and.aspect <-
    function(prepanel.default,
             prepanel = NULL,
             have.xlim = FALSE, xlim = NULL,
             have.ylim = FALSE, ylim = NULL,
             x.relation, y.relation,
             panel.args.common = list(),
             panel.args = list(),
             aspect,
             banking = lattice.getOption("banking"),
             npackets = length(panel.args),
             x.axs = "r", y.axs = "r",
             ...)  ## extra arguments for prepanel (for qqmathline)
{
    prepanel.default.function <- getFunctionOrName(prepanel.default)
    prepanel <- getFunctionOrName(prepanel)
    if (npackets<1) stop("need at least one panel")
    x.limits <- vector("list", npackets)
    y.limits <- vector("list", npackets)
    x.used.at <- vector("list", npackets)
    y.used.at <- vector("list", npackets)
    x.num.limit <- vector("list", npackets)
    y.num.limit <- vector("list", npackets)
    dxdy <- vector("list", npackets)

    for (count in seq_len(npackets))
    {
        if (is.list(panel.args[[count]]))
        {
            pargs <- c(panel.args.common, panel.args[[count]], list(...))
            tem <- do.call("prepanel.default.function", pargs)
            if (is.function(prepanel)) ## results will 'overwrite' defaults
            {
                prenames <- names(formals(prepanel))
                if (!("..." %in% prenames)) pargs <- pargs[intersect(names(pargs), prenames)]
                pretem <- do.call("prepanel", pargs)
                ## prepanel() over-rides defaults except NAs - e.g. ylim = c(0, NA)
                if (!is.null(pretem$xlim) && !is.character(pretem$xlim))
                    if (any(isna <- is.na(pretem$xlim)))
                        pretem$xlim[isna] <- tem$xlim[isna]
                if (!is.null(pretem$ylim) && !is.character(pretem$ylim))
                    if (any(isna <- is.na(pretem$ylim)))
                        pretem$ylim[isna] <- tem$ylim[isna]
                tem <- updateList(tem, pretem)
                ## tem[names(pretem)] <- pretem
            }
            x.limits[[count]] <- tem$xlim
            y.limits[[count]] <- tem$ylim
            x.used.at[[count]] <- if (is.null(tem$xat)) NA else tem$xat
            y.used.at[[count]] <- if (is.null(tem$yat)) NA else tem$yat
            x.num.limit[[count]] <- if (is.null(tem$xat)) NA else range(tem$xat)
            y.num.limit[[count]] <- if (is.null(tem$yat)) NA else range(tem$yat)
            dxdy[[count]] <- list(dx = tem$dx, dy = tem$dy)
        }
        else  ## this happens for empty panels
        {
            x.limits[[count]] <- c(NA_real_, NA_real_)
            y.limits[[count]] <- c(NA_real_, NA_real_)
            x.used.at[[count]] <- NA_real_
            y.used.at[[count]] <- NA_real_
            x.num.limit[[count]] <- NA_real_
            y.num.limit[[count]] <- NA_real_
            dxdy[[count]] <- list(dx = NA_real_, dy = NA_real_)
        }
    }

    ## Some explanation might be helpful here. The for loop above
    ## creates a list of xlims/ylims. Each of these might be either
    ## numeric (when x/y is numeric, shingle or POSIXt etc), or levels
    ## of a factor (that's how prepanel.default.functions are set
    ## up). However, at this point, all x.limits[[i]] must be of the
    ## same type. Returned limits must be in accordance with this
    ## type. The only exception is when relation = "free", in which
    ## case they may be different. This could happen if [xy]lim or
    ## limits is supplied as a list in the high level function.

    x.limits <-
        limitsFromLimitlist(have.lim = have.xlim,
                            lim = xlim,
                            relation = x.relation,
                            limitlist = x.limits,
                            used.at = x.used.at,
                            numlimitlist = x.num.limit,
                            axs = x.axs,
                            npackets = npackets)
    y.limits <-
        limitsFromLimitlist(have.lim = have.ylim,
                            lim = ylim,
                            relation = y.relation,
                            limitlist = y.limits,
                            used.at = y.used.at,
                            numlimitlist = y.num.limit,
                            axs = y.axs,
                            npackets = npackets)

    if (is.character(aspect))
    {
        if (aspect == "xy")
        {
            aspect <-
                median(sapply(dxdy, banking) *
                       y.limits$slicelen /
                       x.limits$slicelen, 
                       na.rm = TRUE)
### old aspect calculation
##             aspect <- median(unlist(lapply(dxdy, banking)),
##                              na.rm = TRUE) * y.limits$slicelen /
##                                  x.limits$slicelen
##             if (y.relation == "free" || x.relation == "free")
##                 warning("'aspect=xy' when 'relation=free' is not sensible")
        }
        else if (aspect == "iso")
        {
            aspect <-
                median(y.limits$slicelen / x.limits$slicelen,
                       na.rm = TRUE)
            if (y.relation == "free" || x.relation == "free")
                warning("'aspect=\"iso\"' approximate since 'relation=\"free\"'")
        }
        else aspect <- 1
    }
    list(x.limits = x.limits$limits,
         y.limits = y.limits$limits,
         x.used.at = x.limits$used.at,
         y.used.at = y.limits$used.at,
         x.num.limit = x.limits$numlimitlist,
         y.num.limit = y.limits$numlimitlist,
         aspect.ratio = aspect,
         prepanel.default = prepanel.default,
         prepanel = prepanel)
}

