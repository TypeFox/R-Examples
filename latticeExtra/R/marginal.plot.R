##
## Copyright (c) 2007 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer

is.categorical <- function (x)
{
    is.factor(x) || is.shingle(x) || is.character(x) || is.logical(x)
}

marginal.plot <-
    function(x,
             data = NULL,
             groups = NULL,
             reorder = !is.table(x),
             plot.points = FALSE,
             ref = TRUE, cut = 0,
             origin = 0, #ylim = c(0, NA), this only supported in R >= 2.11
             xlab = NULL, ylab = NULL,
             type = c("p", if (is.null(groups)) "h"),
             ...,
             subset = TRUE,
             as.table = TRUE,
             subscripts = TRUE,
             default.scales = list(
               relation = "free",
               abbreviate = TRUE, minlength = 5,
               rot = 30, cex = 0.75, tick.number = 3,
               y = list(draw = FALSE)),
             layout = NULL,
             lattice.options = list(
               layout.heights = list(
                 axis.xlab.padding = list(x = 0),
                 xlab.key.padding = list(x = 0))))
{
    if (is.table(data))
        data <- as.data.frame(data)
    ## assume first term of formula is the data object; ignore rest
    if (inherits(x, "formula"))
        x <- eval(x[[2]], data, environment(x))
    ## x must be either a data.frame or a table
    if (!is.data.frame(x) && !is.table(x))
        x <- as.data.frame(x)
    ## groups and subset are subject to non-standard evaluation:
    groups <- eval(substitute(groups), data, parent.frame())
    ## note unusual cases e.g.
    ## evalq(marginal.plot(dat, subset = complete.cases(dat)), myEnv)
    subset <- eval(substitute(subset), data, parent.frame())
    ## apply subset
    if ((length(subset) > 0) && !isTRUE(subset)) {
        x <- x[subset,]
        if (!is.null(groups))
            groups <- groups[subset]
    }
    ## divide into categoricals and numerics
    if (is.table(x)) {
        iscat <- TRUE
    } else {
        iscat <- sapply(x, is.categorical)
    }
    ## reorder factor levels
    if (reorder) {
        if (is.table(x)) {
            x <- reorderTableByFreq(x)
        } else {
            for (nm in names(x)[iscat]) {
                val <- x[[nm]]
                if (is.character(val))
                    x[[nm]] <- factor(val)
                if (!is.ordered(val) &&
                    !is.shingle(val) &&
                    nlevels(val) > 1)
                {
                    x[[nm]] <- reorder(val, val, function(z) -length(z))
                }
            }
        }
    }
    if (any(iscat)) {
        ## handle categorical variables
        ## make a list of dotplot trellis objects
        if (is.table(x)) {
            margins <- seq(length = length(dim(x)))
            names(margins) <- names(dimnames(x))
        } else {
            margins <- which(iscat)
            names(margins) <- colnames(x)[iscat]
        }
        dotobjs <-
            lapply(margins,
                   function(i)
               {
                   if (is.table(x)) {
                       nm <- names(dimnames(x))[i]
                       nm <- deparse(as.symbol(nm), backtick = TRUE)
                       form <- paste("Freq ~", nm)
                       if (!is.null(groups))
                           form <- paste(form, "+ groups")
                       tab <- xtabs(as.formula(form), x)
                   } else {
                       if (!is.null(groups)) {
                           tab <- table(Value = x[[i]], groups = groups)
                       } else {
                           tab <- table(Value = x[[i]])
                       }
                   }
                   dotplot(tab, horizontal = FALSE,
                           groups = !is.null(groups),
                           subscripts = TRUE,
                           ...,
                           type = type,
                           origin = origin, #ylim = ylim,
                           as.table = as.table,
                           default.scales = default.scales,
                           lattice.options = lattice.options,
                           xlab = xlab, ylab = ylab)
               })
        ## merge the list of trellis objects into one
        catobj <- do.call("c", c(dotobjs, merge.legends = FALSE))
        catobj$layout <- layout
        catobj$call <- match.call()
    }
    if (any(!iscat)) {
        ## handle numeric variables
        ## construct formula with all numeric variables
        nms <- names(x)[!iscat]
        symbolStr <- function(nm)
            deparse(as.symbol(nm), backtick = TRUE)
        nms <- sapply(nms, symbolStr)
        numform <- paste("~", paste(nms, collapse = " + "))
        numobj <-
            densityplot(as.formula(numform), x, outer = TRUE,
                        subscripts = TRUE,
                        groups = groups,
                        ...,
                        plot.points = plot.points,
                        ref = ref, cut = cut, #ylim = ylim,
                        as.table = as.table,
                        default.scales = default.scales,
                        lattice.options = lattice.options,
                        xlab = xlab, ylab = ylab)
        ## set strip name if only one panel
        if (prod(dim(numobj)) == 1)
            rownames(numobj) <- names(x)[!iscat]
        numobj$call <- match.call()
        numobj$layout <- layout
    }
    if (all(iscat)) {
        obj <- catobj
    } else if (all(!iscat)) {
        obj <- numobj
    } else {
        ## if there are both categoricals and numerics,
        ## merge the trellis objects; keep original var order
        reIndex <- order(c(which(iscat), which(!iscat)))
        obj <- update(c(catobj, numobj, merge.legends = FALSE),
                      index.cond = list(reIndex), layout = layout)
        ## force strips when only one panel in each object
        if (identical(obj$strip, FALSE))
            obj$strip <- "strip.default"
    }
    obj$call <- sys.call(sys.parent())
    obj
}

reorderTableByFreq <- function(x)
{
    stopifnot(is.table(x))
    df <- as.data.frame(x)
    i <- which(names(df) == "Freq")
    df[-i] <- lapply(df[-i], reorder, - df$Freq)
    xtabs(Freq ~ ., df)
}

