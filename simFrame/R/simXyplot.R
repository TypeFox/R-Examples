# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

setMethod("simXyplot", 
    signature(x = "SimResults"),
    function(x, true = NULL, epsilon, NArate, select, 
            cond = c("Epsilon", "NArate"), 
            average = c("mean", "median"), ...) {
        
        # initializations
        values <- getValues(x)
        if(nrow(values) == 0) stop("slot 'values' does not contain any rows")
        if(ncol(values) == 0) stop("slot 'values' does not contain any columns")
        eps <- getEpsilon(x)
        neps <- length(eps)
        haveEpsilon <- neps > 0
        supEpsilon <- !missing(epsilon)  # contamination level supplied
#        NAr <- getNArate(x)
#        nNAr <- getLength(NAr)
        NAr <- convertNArate(getNArate(x))
        nNAr <- length(NAr)
        haveNArate <- isTRUE(nNAr > 0)
        supNArate <- !missing(NArate)  # missing value rate supplied
        design <- getDesign(x)  # variables used to split simulations
        
        # check contamination levels and missing value rates
#        if(neps <= 1 && isTRUE(nNAr <= 1)) {
#            stop("this plot is only meaningful with either varying ", 
#                "contamination levels or varying missing value rates")
#        } else if(haveEpsilon && haveNArate) {
#            if(neps > 1 && isTRUE(nNAr == 1)) {
#                if(supEpsilon) warning("'epsilon' is ignored")
#                if(supNArate) warning("'NArate' is ignored")
#                xnam <- "Epsilon"
#            } else if(neps == 1 && isTRUE(nNAr > 1)) {
#                if(supEpsilon) warning("'epsilon' is ignored")
#                if(supNArate) warning("'NArate' is ignored")
#                xnam <- "NArate"
#            } else {
#                if(supEpsilon && supNArate) {
#                    stop("only one of 'epsilon' or 'NArate' may be supplied")
#                } else if(supEpsilon) {
#                    if(!isTRUE(length(epsilon) == 1)) {
#                        stop("'epsilon' must specify ", 
#                        "exactly one contamination level")
#                    }
#                    values <- values[values$Epsilon == epsilon, , drop=FALSE]
#                    if(nrow(values) == 0) {
#                        stop("the contamination level specified by ",
#                            "'epsilon' has not been used in the simulation")
#                    }
#                    xnam <- "NArate"
#                } else if(supNArate) {
#                    if(!isTRUE(length(NArate) == 1)) {
#                        stop("'NArate' must specify ", 
#                            "exactly one missing value rate")
#                    }
#                    values <- values[values$NArate == NArate, , drop=FALSE]
#                    if(nrow(values) == 0) {
#                        stop("the missing value rate specified by ",
#                            "'NArate' has not been used in the simulation")
#                    }
#                    xnam <- "Epsilon"
#                } else {
#                    stop("contamination levels and missing ", 
#                        "value rates are both varying")
#                }
#            }
#        } else if(haveEpsilon) {
#            if(supEpsilon) warning("'epsilon' is ignored")
#            xnam <- "Epsilon"
#        } else if(haveNArate) {
#            if(supNArate) warning("'NArate' is ignored")
#            xnam <- "NArate"
#        } else stop("unexpected problem with 'x'")  # just to be safe
        if(haveEpsilon && supEpsilon) {
            if(neps == 1) warning("'epsilon' is ignored")
            else {
                eps <- intersect(eps, epsilon)
                neps <- length(eps)
                if(neps == 0) {
                    stop("the contamination levels specified by ",
                        "'epsilon' have not been used in the simulation")
                } else values <- values[values$Epsilon %in% eps, , drop=FALSE]
            }
        }
        if(haveNArate && supNArate) {
            if(nNAr == 1) warning("'NArate' is ignored")
            else {
                NAr <- intersect(NAr, NArate)
                nNAr <- length(NAr)
                if(nNAr == 0) {
                    stop("the missing value rates specified by ",
                        "'NArate' have not been used in the simulation")
                } else values <- values[values$NArate %in% NAr, , drop=FALSE]
            }
        } 
        if(neps <= 1 && isTRUE(nNAr <= 1)) {
            stop("this plot is only meaningful with either varying ", 
                "contamination levels or varying missing value rates")
        } else if(haveEpsilon && haveNArate) {
            if(neps > 1 && nNAr == 1) xnam <- "Epsilon"
            else if(neps == 1 && nNAr > 1) xnam <- "NArate"
            else {
                if(missing(cond)) {
                    if(neps < nNAr) cond <- "Epsilon"
                    else cond <- "NArate"
                } else cond <- match.arg(cond)
                xnam <- setdiff(c("Epsilon", "NArate"), cond)
                design <- c(design, cond)
                f <- as.factor(paste(cond, "=", values[, cond]))
#                levels(f) <- rev(levels(f))
                values[, cond] <- f
                if(!is.null(true)) true <- rep(true, length(levels(f)))
            }
        } else if(haveEpsilon) xnam <- "Epsilon"
        else if(haveNArate) xnam <- "NArate"
        else stop("unexpected problem with 'x'")  # just to be safe
        
        # check specified columns
        if(missing(select)) select <- getColnames(x)
        else {
            if(!is.character(select)) {
                stop("'select' must be a character vector")
            }
            if(!all(select %in% getColnames(x))) stop("undefined columns selected")
        }
        
        # check method for computing average
        average <- match.arg(average)
        
        # if missing value rates are plotted on x-axis and NArate is a matrix, 
        # the display on the x-axis should be more of a categorical nature 
        # (corresponding to the rows of NArate)
#        at <- if(haveNArate && is(getNArate(x), "matrix")) 1:nNAr else NULL
        at <- if(xnam == "NArate" && is(getNArate(x), "matrix")) NAr else NULL

        # call internal function
        internalSimXyplot(values, xnam, design, select, 
            at=at, true=true, average=average, ...)
    })


## internal function
internalSimXyplot <- function(values, xnam = c("Epsilon", "NArate"), 
        design, names, at = NULL, true = NULL, average, auto.key=TRUE, 
        scales=list(), type = "l", ylab = NULL, ..., 
        # the following arguments are defined so that they aren't supplied twice
        x, data, panel, prepanel, groups) {
    # prepare legend
    if(isTRUE(auto.key)) auto.key <- list(points=FALSE, lines=TRUE)
    else if(is.list(auto.key)) {
        if(is.null(auto.key$points)) auto.key$points <- FALSE
        if(is.null(auto.key$lines)) auto.key$lines <- TRUE
    }
    # prepare data and formula for lattice graphics
    xnam <- match.arg(xnam)
    if(xnam == "NArate") {
        if(!is.null(at)) {
            if(is.null(scales$x$at)) scales$x$at <- at
            if(is.null(scales$x$tck)) scales$x$tck <- 0
        }
    }
    cond <- c(xnam, design)
    if(average == "median") f <- function(x) median(x, na.rm=TRUE)
    else f <- function(x) mean(x, na.rm=TRUE)
    tmp <- aggregate(values[, names, drop=FALSE], values[, cond, drop=FALSE], f)
    values <- getLatticeData(tmp, cond, names)
    form <- getFormula(".Value", xnam, design)
    # call lattice function
    # this produces a 'NOTE' during 'R CMD check':
#    xyplot(form, data=values, auto.key=auto.key, panel=panelSimXyplot, 
#        prepanel=prepanelSimXyplot, scales=scales, true=true, groups=.Name, 
#        type=type, ylab=ylab, ...)
    command <- paste("xyplot(form, data=values, auto.key=auto.key,", 
        "panel=panelSimXyplot, prepanel=prepanelSimXyplot, scales=scales,", 
        "true=true, groups=.Name, type=type, ylab=ylab, ...)")
    eval(parse(text=command))
}

## panel function
panelSimXyplot <- function(x, y, true = NULL, ...) {
    if(length(true) > 0) {
        # if only one value is supplied for true, it is used for 
        # every panel, otherwise we select the one corresponding
        # to the panel
        if(length(true) > 1) {
            i <- packet.number()  # identifies current packet
            true <- true[i]
        }
        panel.refline(h=true, ...)
    }
    panel.xyplot(x, y, ...)
}

## prepanel function
prepanelSimXyplot <- function(x, y, true = NULL, ...) {
    if(length(true) == 0) list()
    else {
# FIXME: default axis limits if relation is "sliced" or "free"
#        'packet.number' only works for panel function
#        can packet number be accessed from the prepanel function?
#        # if only one value is supplied for true, it is used for 
#        # every panel, otherwise we select the one corresponding
#        # to the panel
#        if(length(true) > 1) {
#            i <- packet.number()  # identifies current packet
#            true <- true[i]
#        }
        lim <- range(c(y, true), finite=TRUE)
        list(ylim=lim)
    }
}
