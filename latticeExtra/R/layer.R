##
## Copyright (c) 2008 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer

as.layer <- function(x, ...)
    UseMethod("as.layer")

as.layer.layer <- function(x, ...)
    x

layer <-
    function(..., data = NULL,
             magicdots = TRUE, exclude = NULL,
             packets = NULL, 
             rows = NULL, columns = NULL,
             groups = NULL,
             style = NULL, force = FALSE,
             theme = if (force) trellis.par.get() else NULL,
             under = FALSE, superpose = FALSE)
{
    ## set layer to quoted expressions in `...`
    foo <- eval(substitute(expression(...)))
    if (magicdots) {
        ## The dots `...` are magic:
        ## pass on only those arguments not named in each call
        foo <- as.expression(lapply(foo, magicDots, exclude = exclude))
    }
    mostattributes(foo) <-
        list(data = data,
             under = under,
             packets = packets,
             rows = rows,
             columns = columns,
             groups = groups,
             superpose = superpose,
             style = style,
             theme = theme)
    lay <- list(foo)
    class(lay) <- c("layer", "trellis")
    lay
}

## convert a call containing `...` to only pass on arguments
## not named in the call
magicDots <- function(ocall, exclude = NULL, assume.xy = TRUE)
{
    if (!is.call(ocall)) stop("arguments to layer() should be calls")
    ## call recursively with any calls inside this one
    for (i in seq_along(ocall)[-1]) {
        thisArg <- ocall[[i]]
        if (missing(thisArg)) ## eg x[,1]
            next
        if (is.call(thisArg)) {
            ## skip function definitions
            if (identical(thisArg[[1]], as.symbol("function")))
                next
            ocall[[i]] <- Recall(thisArg, exclude = exclude, assume.xy = assume.xy)
        }
    }
    Args <- as.list(ocall)[-1]
    ## nothing to do if there are no dots in the call
    idots <- sapply(Args, identical, as.symbol("..."))
    if (!any(idots))
        return(ocall)
    Args <- Args[!idots]
    ## nothing to do if there are only dots in the call (unless exclude)
    if ((length(Args) == 0) && (length(exclude) == 0))
        return(ocall)
    ## assume first argument is 'x' if is un-named, and second 'y'
    if (assume.xy && (length(Args) > 0)) {
        if (is.null(names(Args)))
            names(Args) <- rep("", length = length(Args))
        if (identical(names(Args)[1], ""))
            names(Args)[1] <- "x"
        if (identical(names(Args)[2], ""))
            names(Args)[2] <- "y"
    }
    if (length(exclude) == 0) {
        ## simple case
        mcall <-
            substitute(do.call(FUN,
                               modifyList(list(...), Args)),
                       list(FUN = ocall[[1]], Args = Args))
    } else {
        ## exclude named arguments from dots
        mcall <-
            substitute(do.call(FUN,
                               modifyList(list(...)[!(names(list(...)) %in% exclude)],
                                          Args)),
                       list(FUN = ocall[[1]], Args = Args, exclude = exclude))
    }
    mcall
}

layer_ <- function(...)
{
    ccall <- match.call()
    ccall$under <- TRUE
    ccall[[1]] <- quote(layer)
    eval.parent(ccall)
}

glayer <- function(...)
{
    ccall <- match.call()
    ccall$superpose <- TRUE
    ccall[[1]] <- quote(layer)
    eval.parent(ccall)
}

glayer_ <- function(...)
{
    ccall <- match.call()
    ccall$superpose <- TRUE
    ccall$under <- TRUE
    ccall[[1]] <- quote(layer)
    eval.parent(ccall)
}

## to avoid print.trellis
print.layer <- function(x, ...) print.default(x, ...)

## to avoid [.trellis and to keep the class attribute
"[.layer" <- function (x, i, ...)
    structure(unclass(x)[i], class = class(x))

"+.trellis" <- function(object, lay)
{
    ocall <- sys.call(sys.parent()); ocall[[1]] <- quote(`+`)
    if (missing(object) || missing(lay)) stop("Only one argument supplied to binary operator + which requires two.")
    stopifnot(inherits(object, "trellis"))
    lay <- as.layer(lay)
    if (inherits(object, "layer")) {
        ## just concatenate lists
        return(structure(c(unclass(object), unclass(lay)),
                         class = c("layer", "trellis")))
    }
    panel <- if ("panel" %in% names(object$panel.args.common))
        object$panel.args.common$panel
    else object$panel
    panel <- if (is.function(panel)) panel
    else if (is.character(panel)) {
        ## could be just get(panel), but for flattenPanel:
        ## do not expand original panel function eg panel.xyplot(...)
        tmp <- function(...) NA
        body(tmp) <- call(panel, quote(...))
        environment(tmp) <- globalenv()
        tmp
    } else eval(panel)
    ## a flag to indicate this panel function has layers
    ## (used by flattenPanel and undoLayer)
    .is.a.layer <- TRUE
    newpanel <- function(...) {
        .UNDER <- unlist(lapply(lay, attr, "under"))
        ## underlaying items only
        drawLayer(lay[.UNDER], list(...))
        ## original panel function:
        panel(...)
        ## overlaying items only
        drawLayer(lay[.UNDER == FALSE], list(...))
    }
    if ("panel" %in% names(object$panel.args.common))
        object$panel.args.common$panel <- newpanel
    else object$panel <- newpanel
    ## need this to allow further calls to update() to insert arguments:
    object$call <- call("update", ocall)
    object
}

drawLayer <- function(lay, panelArgs = trellis.panelArgs())
{
    lay <- as.layer(lay)
    .UNDER <- unlist(lapply(lay, attr, "under"))
    ## underlayers, in reverse order
    for (.ITEM in rev(lay[.UNDER]))
        drawLayerItem(.ITEM, panelArgs)
    ## overlayers
    for (.ITEM in lay[.UNDER == FALSE])
        drawLayerItem(.ITEM, panelArgs)
    invisible()
}

drawLayerItem <- function(layer.item, panelArgs)
{
    stopifnot(is.expression(layer.item))
    ## check that any restrictions on packets/rows/columns are met
    matchesok <- function(spec, value) {
        if (is.null(spec)) return(TRUE)
        if (is.numeric(spec) && all(spec <= 0))
            ## negative indexes exclude items
            return(value %in% -spec == FALSE)
        else
            return(value %in% spec)
    }
    matchesallok <-
        with(list(a = attributes(layer.item)),
             matchesok(a$packets, packet.number()) &&
             matchesok(a$rows, current.row()) &&
             matchesok(a$columns, current.column()))
    if (!matchesallok) return()
    ## set given theme for duration of this function
    if (!is.null(attr(layer.item, "theme"))) {
        .TRELLISPAR <- trellis.par.get()
        trellis.par.set(attr(layer.item, "theme"))
        on.exit(trellis.par.set(.TRELLISPAR))
    }
    ## define a layer drawing function, which may be per group
    drawLayerItemPerGroup <- function(...)
    {
        ## Note: layer.item is found in this function's environment
        dots <- list(...)
        ## restrict to specified group numbers
        groupok <- (matchesok(attr(layer.item, "groups"), dots$group.number) ||
                    matchesok(attr(layer.item, "groups"), as.character(dots$group.value)))
        if (!groupok)
            return()
        if (!is.null(attr(layer.item, "style"))) {
            ## extract plot style attributes from given index into superpose.*
            .TRELLISPAR <- trellis.par.get()
            local({
                i <- attr(layer.item, "style")
                line <- Rows(trellis.par.get("superpose.line"), i)
                symbol <- Rows(trellis.par.get("superpose.symbol"), i)
                polygon <- Rows(trellis.par.get("superpose.polygon"), i)
                trellis.par.set(plot.line = line,
                                superpose.line = line,
                                add.line = line,
                                add.text = line,
                                plot.symbol = symbol,
                                superpose.symbol = symbol,
                                plot.polygon = polygon,
                                superpose.polygon = polygon,
                                axis.text = line,
                                axis.line = line
                                )
            })
            on.exit(trellis.par.set(.TRELLISPAR))
        }
        with(dots,
             eval(layer.item, attr(layer.item, "data"),
                  environment()))
    }
    ## call panel.superpose for group layers
    if (isTRUE(attr(layer.item, "superpose"))) {
        do.call("panel.superpose",
                modifyList(panelArgs,
                  list(panel.groups = drawLayerItemPerGroup)))
    } else {
        do.call("drawLayerItemPerGroup", panelArgs)
    }
}

flattenPanel <- function(object)
{
    flattenFun <- function(fun)
    {
        env <- environment(fun)
        ## check if this panel function is simple or has layers
        if (is.null(env) ||
            !exists(".is.a.layer", env, inherits = FALSE))
            return(as.expression(body(fun)))
        ## merge: under layers, existing panel, over layers
        .UNDER <- sapply(env$lay, attr, "under")
        c(do.call("c", rev(env$lay[.UNDER])),
          flattenFun(env$panel),
          do.call("c", env$lay[.UNDER == FALSE]))
    }
    flat <- flattenFun(object$panel)
    ## wrap in braces, as in a function body
    as.call(c(quote(`{`), flat))
}

## not exported -- I do not think this is really useful
undoLayer <- function(x)
{
    stopifnot(is.function(x$panel))
    env <- environment(x$panel)
    if (!exists(".is.a.layer", env, inherits=FALSE))
        stop("does not look like a layer")
    update(x, panel=env$panel)
}

