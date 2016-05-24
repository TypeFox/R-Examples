# playwith: interactive plots in R using GTK+
##
## Copyright (c) 2007 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer

getDataArg <- function(playState = playDevCur(), eval = TRUE)
{
    tmp.data <- NULL
    if (!is.null(playState$data.points)) {
        ## data.points were supplied
        tmp.data <- if (eval) playState$data.points
        else quote(playDevCur()$data.points) ## unusual
    } else {
        ## data.points missing; guess from call
        mainCall <- mainCall(playState)
        if (length(mainCall) > 1) {
            ## check for named "data" argument
            tmp.data <- callArg(playState, "data", eval = eval)
            if (is.null(tmp.data)) {
                ## look at first argument
                tmp.x <- callArg(playState, 1, eval = TRUE)
                if (inherits(tmp.x, "formula")) {
                    ## if 1st arg is formula, 2nd is data if un-named (by convention)
                    if (is.null(tmp.data) &&
                        (length(mainCall) >= 3) &&
                        (is.null(names(mainCall)) ||
                         identical(names(mainCall)[[3]], "")))
                    {
                        tmp.data <- callArg(playState, 2, eval = eval)
                    }
                } else if (!is.vector(tmp.x)) {
                    ## take first argument as data object, but need to be careful:
                    ## can we coerce to data.frame? (e.g. not for "density")
                    ## need to coerce e.g. sp::SpatialPoints to data.frame
                    tmp.x.df <- try(as.data.frame(tmp.x), silent = TRUE)
                    if (!inherits(tmp.x.df, "try-error")) {
                        if (eval) {
                            tmp.data <- tmp.x.df
                        } else {
                            tmp.data <- callArg(playState, 1, eval = FALSE)
                            tmp.data <- call("as.data.frame", tmp.data)
                        }
                    }
                }
            }
        }
        if (is.null(tmp.data)) {
            ## objects may also come from a with() block
            if (identical(playState$call[[1]], as.symbol("with"))) {
                tmp.data <- playState$call[[2]]
                if (eval) tmp.data <- eval(tmp.data, playState$env)
            }
        }
    }
    tmp.data
}

guessLabels <- function(playState)
{
    playState$tmp$identify.ok <- FALSE
    if (isBasicDeviceMode(playState)) return()
    hasArgs <- playState$accepts.arguments
    isLatt <- playState$is.lattice
    isBase <- !isLatt && is.null(playState$viewport)
    isBaseMulti <- isBase && any(par("mfrow") > 1)
    ## detect known plots that this will not work with
    if (is.null(playState$data.points)) {
        if (!hasArgs) return()
        if (isBaseMulti) return()
        ## lattice package:
        if (playState$callName %in%
            c("histogram", "densityplot", "barchart")) return()
        ## latticeExtra package:
        if (playState$callName %in% "marginal.plot") return()
        ## graphics package:
        if (playState$callName %in%
            c("hist", "barplot", "spineplot", "mosaic",
              "assocplot", "fourfoldplot",
              "coplot", "persp", "pie")) return()
    }
    ok.case.names <- function(object) {
        tmp <- case.names(object)
        if (is.null(tmp))
            return(seq_len(NROW(object)))
        tmp
    }
    ## try to guess labels if they were not given
    labels <- playState$.args$labels
    if (is.null(labels)) {
        tmp.data <- getDataArg(playState)
        if (!is.null(tmp.data) &&
            !inherits(tmp.data, "list") &&
            !is.environment(tmp.data))
        {
            ## data arg, probably a data.frame
            labels <- ok.case.names(tmp.data)
        } else {
            ## no useful data arg; take arg 1 instead
            tmp.x <- callArg(playState, 1)
            if (inherits(tmp.x, "formula")) {
                ## get left-most term in RHS of formula
                xObj <- if (length(tmp.x) == 2)
                    tmp.x[[2]] else tmp.x[[3]]
                while (is.call(xObj) && toString(xObj[[1]]) %in%
                       c("|", "*", "+"))
                    xObj <- xObj[[2]]
                xObj <- eval(xObj, tmp.data, playState$env)
                labels <- ok.case.names(xObj)
            } else {
                ## first arg is an object, not a formula
                labels <- ok.case.names(tmp.x)
            }
        }
    }
    playState$labels <- labels
    playState$tmp$identify.ok <- TRUE
    invisible()
}

xyCoords <- function(playState = playDevCur(), space="plot")
{
    foo <- xyData(playState, space=space)
    ## factor data can't be stored as a matrix;
    ## if either is a matrix, convert
    ## both x and y to same matrix form.
    nvar <- max(NCOL(foo$x), NCOL(foo$y))
    ## if 'subscripts' exists, that can also
    ## define the number of cases
    if (length(foo$subscripts) > 1) {
        nvar <- min(nvar, length(foo$subscripts))
    }
    if (!is.numeric(foo$x)) {
        foo$x <- as.numeric(foo$x)
        if (nvar > 1)
            foo$x <- matrix(foo$x, ncol = nvar)
    }
    if (!is.numeric(foo$y)) {
        foo$y <- as.numeric(foo$y)
        if (nvar > 1)
            foo$y <- matrix(foo$y, ncol = nvar)
    }
    foo
}

xyData <- function(playState = playDevCur(), space="plot")
{
    if (!is.null(playState$data.points)) {
        return(xy.coords_with_class(playState$data.points))
    }
    ## pass call name with a self-titled class to allow user-defined
    ## plotCoords definitions
    name <- playState$callName
    class(name) <- name
    call <- mainCall(playState)
    envir <- playState$env

    ## lattice -- get panel data from stored trellis object
    if (playState$is.lattice) {
        if (space == "page") {
            ## TODO: this only works with standard coords (xyplot etc)
            ## data from all packets
            tmp <- try(do.call(rbind,
                               lapply(playState$trellis$panel.args, as.data.frame)
                               ), silent=TRUE)
            if (inherits(tmp, "try-error"))
                return(NULL)
            return(tmp)
        }
        if (space == "plot") {
            space <- packet.number()
            if (length(space) == 0) {
                packets <- playState$tmp$currentLayout
                if (sum(packets > 0) > 1) stop("space not well specified")
                space <- packets[packets > 0][1]
            }
            space <- paste("packet", space)
        }
        packet <- as.numeric(sub("packet ", "", space))
        if (is.na(packet) ||
            (packet > length(playState$trellis$panel.args)))
            return(NULL)
        panel.args <- trellis.panelArgs(playState$trellis, packet)
        ## look for plotCoords method and pass panel.args
        object <- if (length(call) > 1)
            callArg(playState, 1)
        return(plotCoords(name, object = object, call = call,
                          envir = envir, panel.args = panel.args,
                          packet = packet))
    } else {
        ## non-lattice plot:
        ## look for plotCoords method and pass data arg
        data <- getDataArg(playState)
        object <- if (length(call) > 1)
            callArg(playState, 1, data = data)
        plotCoords(name, object = object, call = call,
                   envir = envir, data = data)
    }
}


## aka playUnLog
spaceCoordsToDataCoords <- function(playState, xy)
{
    if (!is.null(xy$x)) xy$x <-
        spaceCoordsToDataCoordsXY(playState, xy$x, x.or.y="x")
    if (!is.null(xy$y)) xy$y <-
        spaceCoordsToDataCoordsXY(playState, xy$y, x.or.y="y")
    xy
}

spaceCoordsToDataCoordsXY <- function(playState, x, x.or.y=c("x", "y"))
{
    x.or.y <- match.arg(x.or.y)
    logBase <- playLogBase(playState, x.or.y)
    if (is.na(logBase)) return(x)
    logBase ^ x
}

## aka playReLog
dataCoordsToSpaceCoords <- function(playState, xy)
{
    if (!is.null(xy$x)) xy$x <-
        dataCoordsToSpaceCoordsXY(playState, xy$x, x.or.y="x")
    if (!is.null(xy$y)) xy$y <-
        dataCoordsToSpaceCoordsXY(playState, xy$y, x.or.y="y")
    xy
}

dataCoordsToSpaceCoordsXY <- function(playState, x, x.or.y=c("x", "y"))
{
    x.or.y <- match.arg(x.or.y)
    logBase <- playLogBase(playState, x.or.y)
    if (is.na(logBase)) return(x)
    log(x, base=logBase)
}

playLogBase <- function(playState, x.or.y=c("x", "y"))
{
    x.or.y <- match.arg(x.or.y)
    if (playState$is.lattice) {
        scalesObj <- playState$trellis[[paste(x.or.y, "scales", sep=".")]]
        x <- scalesObj$log
        if (identical(x, FALSE)) return(NA)
        if (isTRUE(x)) return(10)
        if (identical(x, "e")) return(exp(1))
        return(x)
    } else if (playState$is.ggplot) {
        logArg <- callArg(playState, "log")
        if (!is.null(logArg) &&
            (x.or.y %in% strsplit(logArg, split="")[[1]]))
            return(10)
    } else {
        ## traditional graphics plot
        if (par(paste(x.or.y, "log", sep=""))) return(10)
    }
    return(NA)
}
