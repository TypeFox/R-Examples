# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

setMethod("simDensityplot", 
    signature(x = "SimResults"),
    function(x, true = NULL, epsilon, NArate, select, ...) {
        
        # initializations
        values <- getValues(x)
        if(nrow(values) == 0) stop("slot 'values' does not contain any rows")
        if(ncol(values) == 0) stop("slot 'values' does not contain any columns")
#        haveEpsilon <- length(getEpsilon(x)) > 0
        eps <- getEpsilon(x)
        neps <- length(eps)
        haveEpsilon <- neps > 0
        supEpsilon <- !missing(epsilon)  # contamination level supplied
#        haveNArate <- isTRUE(getLength(getNArate(x)) > 0)
        NAr <- convertNArate(getNArate(x))
        nNAr <- length(NAr)
        haveNArate <- isTRUE(nNAr > 0)
        supNArate <- !missing(NArate)  # missing value rate supplied
        design <- getDesign(x)  # variables used to split simulations
        
        # check contamination level and missing value rate
#        if(haveEpsilon) {
#            if(supEpsilon) {
#                if(!isTRUE(length(epsilon) == 1)) {
#                    stop("'epsilon' must specify ", 
#                        "exactly one contamination level")
#                }
#                eps <- values$Epsilon
#            } else if(length(getEpsilon(x)) > 1) {
#                stop("this plot is not meaningful ", 
#                    "for varying contamination levels")
#            }
#        }
#        if(haveNArate) {
#            if(supNArate) {
#                if(!isTRUE(length(NArate) == 1)) {
#                    stop("'NArate' must specify ", 
#                        "exactly one missing value rate")
#                }
#                NAr <- values$NArate
#            } else if(isTRUE(getLength(getNArate(x)) > 1)) {
#                stop("this plot is not meaningful ", 
#                    "for varying missing value rates")
#            }
#        }
#        # get indices of selected results
#        if(haveEpsilon && supEpsilon && haveNArate && supNArate) {
#            values <- values[(eps == epsilon) & (NAr == NArate), , drop=FALSE]
#            if(nrow(values) == 0) {
#                stop("the contamination level specified by ",
#                    "'epsilon' or the missing value rate specified by ", 
#                    "'NArate' has not been used in the simulation")
#            }
#        } else if(haveEpsilon && supEpsilon) {
#            values <- values[eps == epsilon, , drop=FALSE]
#            if(nrow(values) == 0) {
#                stop("the contamination level specified by ",
#                    "'epsilon' has not been used in the simulation")
#            }
#        } else if(haveNArate && supNArate) {
#            values <- values[NAr == NArate, , drop=FALSE]
#            if(nrow(values) == 0) {
#                stop("the missing value rate specified by ",
#                    "'NArate' has not been used in the simulation")
#            }
#        }
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
        ## use NA rate first for conditioning in lattic panels
        if(nNAr > 1) {
            design <- c(design, "NArate")
            f <- as.factor(paste("NArate", "=", values[, "NArate"]))
#            levels(f) <- rev(levels(f))
            values[, "NArate"] <- f
            if(!is.null(true)) true <- rep(true, length(levels(f)))
        }
        ## use epsilon as top level for conditioning in lattic panels
        if(neps > 1) {
            design <- c(design, "Epsilon")
            f <- as.factor(paste("Epsilon", "=", values[, "Epsilon"]))
#            levels(f) <- rev(levels(f))
            values[, "Epsilon"] <- f
            if(!is.null(true)) true <- rep(true, length(levels(f)))
        }
        
        # check specified columns
        if(missing(select)) select <- getColnames(x)
        else {
            if(!is.character(select)) {
                stop("'select' must be a character vector")
            }
            if(!all(select %in% getColnames(x))) stop("undefined columns selected")
        }
        
        # call internal function
        internalSimDensityplot(values, design, select, true=true, ...)
    })


## internal function
internalSimDensityplot <- function(values, design, names, 
        true = NULL, auto.key=TRUE, xlab = NULL, ..., 
        # the following arguments are defined so that they aren't supplied twice
        x, data, panel, prepanel, groups) {
    # prepare legend
    if(isTRUE(auto.key)) auto.key <- list(points=TRUE, lines=TRUE)
    else if(is.list(auto.key)) {
        if(is.null(auto.key$points)) auto.key$points <- TRUE
        if(is.null(auto.key$lines)) auto.key$lines <- TRUE
    }
    # prepare data and formula for lattice graphics
    values <- getLatticeData(values, design, names)
    form <- getFormula("", ".Value", design)
    # call lattice function
    # this produces a 'NOTE' during 'R CMD check':
#    densityplot(form, data=values, auto.key=auto.key, 
#        panel=panelSimDensityplot, prepanel=prepanelSimDensityplot, 
#        true=true, groups=.Name, xlab=xlab, ...)
    command <- paste("densityplot(form, data=values, auto.key=auto.key,",  
        "panel=panelSimDensityplot, prepanel=prepanelSimDensityplot,",  
        "true=true, groups=.Name, xlab=xlab, ...)")
    eval(parse(text=command))
}

## panel function
panelSimDensityplot <- function(x, true = NULL, ...) {
    if(length(true) > 0) {
        # if only one value is supplied for true, it is used for 
        # every panel, otherwise we select the one corresponding
        # to the panel
        if(length(true) > 1) {
            i <- packet.number()  # identifies current packet
            true <- true[i]
        }
        panel.refline(v=true, ...)
    }
    panel.densityplot(x, ...)
}

## prepanel function
prepanelSimDensityplot <- function(x, true = NULL, ...) {
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
        lim <- range(c(x, true), finite=TRUE)
        list(xlim=lim)
    }
}
