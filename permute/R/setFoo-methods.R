## Replacement functions for blocks, plots and within, plus strata,
## etc ...
`setNperm<-` <- function(object, value) {
    UseMethod("setNperm<-")
}

`setNperm<-.default` <- function(object, value) {
    stop("No default method for `setNperm`")
}

`setNperm<-.how` <- function(object, value) {
    object[["nperm"]] <- value
    object <- fixupCall(object, "nperm", value)
    object
}

`setMaxperm<-` <- function(object, value) {
    UseMethod("setMaxperm<-")
}

`setMaxperm<-.default` <- function(object, value) {
    stop("No default method for `setMaxperm`")
}

`setMaxperm<-.how` <- function(object, value) {
    object[["maxperm"]] <- value
    object <- fixupCall(object, "maxperm", value)
    object
}

`setMinperm<-` <- function(object, value) {
    UseMethod("setMinperm<-")
}

`setMinperm<-.default` <- function(object, value) {
    stop("No default method for `setMinperm`")
}

`setMinperm<-.how` <- function(object, value) {
    object[["minperm"]] <- value
    object <- fixupCall(object, "minperm", value)
    object
}

`setComplete<-` <- function(object, value) {
    UseMethod("setComplete<-")
}

`setComplete<-.default` <- function(object, value) {
    stop("No default method for `setComplete`")
}

`setComplete<-.how` <- function(object, value) {
    if (!is.null(value))
        value <- rep(as.logical(value), length.out = 1)
    object[["complete"]] <- value
    object <- fixupCall(object, "complete", value)
    object
}

`setAllperms<-` <- function(object, value) {
    UseMethod("setAllperms<-")
}

`setAllperms<-.default` <- function(object, value) {
    stop("No default method for `setAllperms`")
}

`setAllperms<-.how` <- function(object, value) {
    if (!is.null(value))
        value <- as.matrix(value)
    object[["all.perms"]] <- value
    object <- fixupCall(object, "all.perms", value)
    object
}

`setMake<-` <- function(object, value) {
    UseMethod("setMake<-")
}

`setMake<-.default` <- function(object, value) {
    stop("No default method for `setMake`")
}

`setMake<-.how` <- function(object, value) {
    if (!is.null(value))
        value <- rep(as.logical(value), length.out = 1)
    object[["make"]] <- value
    object <- fixupCall(object, "make", value)
    object
}

`setBlocks<-` <- function(object, value) {
    UseMethod("setBlocks<-")
}

`setBlocks<-.default` <- function(object, value) {
    stop("No default method for `setBlocks`")
}

`setBlocks<-.how` <- function(object, value) {
    object[["blocks.name"]] <- deparse(substitute(value))
    if (!is.null(value))
        value <- as.factor(value)
    object["blocks"] <- list(value)
    object <- fixupCall(object, "blocks", value)
    object
}

`setObserved<-` <- function(object, value) {
    UseMethod("setObserved<-")
}

`setObserved<-.default` <- function(object, value) {
    stop("No default method for `setObserved`")
}

`setObserved<-.how` <- function(object, value) {
    if (!is.null(value))
        value <- rep(as.logical(value), length.out = 1)
    object[["observed"]] <- value
    object <- fixupCall(object, "observed", value)
    object
}

## Plots ##############################################################
`setPlots<-` <- function(object, value) {
    UseMethod("setPlots<-")
}

`setPlots<-.default` <- function(object, value) {
    stop("No default method for `setPlots`")
}

`setPlots<-.how` <- function(object, value) {
    stopifnot(inherits(value, "Plots"))
    object[["plots"]] <- value
    object <- fixupCall(object, "plots", getCall(value))
    object
}

## Within ##############################################################
`setWithin<-` <- function(object, value) {
    UseMethod("setWithin<-")
}

`setWithin<-.default` <- function(object, value) {
    stop("No default method for `setWithin`")
}

`setWithin<-.how` <- function(object, value) {
    stopifnot(inherits(value, "Within"))
    object[["within"]] <- value
    object <- fixupCall(object, "within", getCall(value))
    object
}

## Strata #############################################################
`setStrata<-` <- function(object, value) {
    UseMethod("setStrata<-")
}

`setStrata<-.default` <- function(object, value) {
    stop("No default method for `setStrata`")
}

`setStrata<-.how` <- function(object, value) {
    if (!is.null(value)) {
        value <- as.factor(value)
    }
    ## get Plots
    plots <- getPlots(object)
    setStrata(plots) <- value
    setPlots(object) <- plots
    object
}

`setStrata<-.Plots` <- function(object, value) {
    if (!is.null(value))
        value <- as.factor(value)
    object[["strata"]] <- value
    object <- fixupCall(object, "strata", value) # value was getCall(value))
    object
}

## Grid dimensions ####################################################
`setRow<-` <- function(object, value) {
    UseMethod("setRow<-")
}

`setRow<-.default` <- function(object, value) {
    stop("No default method for `setRow`")
}

`setRow<-.how` <- function(object, value) {
    stop("`setRow` can not be used directly on '\"how\"' objects.")
}

`setRow<-.Within` <- function(object, value) {
    value <- as.integer(value)
    object[["nrow"]] <- value
    object <- fixupCall(object, "nrow", value)
    object
}

`setRow<-.Plots` <- function(object, value) {
    value <- as.integer(value)
    object[["nrow"]] <- value
    object <- fixupCall(object, "nrow", value)
    object
}

`setCol<-` <- function(object, value) {
    UseMethod("setCol<-")
}

`setCol<-.default` <- function(object, value) {
    stop("No default method for `setCol`")
}

`setCol<-.how` <- function(object, value) {
    stop("`setCol` can not be used directly on '\"how\"' objects.")
}

`setCol<-.Within` <- function(object, value) {
    value <- as.integer(value)
    object[["ncol"]] <- value
    object <- fixupCall(object, "ncol", value)
    object
}

`setCol<-.Plots` <- function(object, value) {
    value <- as.integer(value)
    object[["ncol"]] <- value
    object <- fixupCall(object, "ncol", value)
    object
}

`setDim<-` <- function(object, value) {
    UseMethod("setDim<-")
}

`setDim<-.default` <- function(object, value) {
    stop("No default method for `setDim`")
}

`setDim<-.how` <- function(object, value) {
    stop("`setDim` can not be used directly on '\"how\"' objects.")
}

`setDim<-.Within` <- function(object, value) {
    value <- as.integer(value)
    stopifnot(all.equal(length(value), 2L))
    setRow(object) <- value[1]
    setCol(object) <- value[2]
    object
}

`setDim<-.Plots` <- function(object, value) {
    value <- as.integer(value)
    stopifnot(all.equal(length(value), 2L))
    setRow(object) <- value[1]
    setCol(object) <- value[2]
    object
}

## setType ############################################################
`setType<-` <- function(object, value) {
    UseMethod("setType<-")
}

`setType<-.default` <- function(object, value) {
    stop("No default method for `setType`")
}

`setType<-.how` <- function(object, value) {
    stop("`setType` can not be used directly on '\"how\"' objects.")
}

`setType<-.Within` <- function(object, value) {
    value <- as.character(value)
    if (!value %in% c("free","series","grid","none"))
        stop("Invalid permutation type")
    value <- rep(value, length.out = 1L)
    object[["type"]] <- value
    object <- fixupCall(object, "type", value)
    object
}

`setType<-.Plots` <- function(object, value) {
    value <- as.character(value)
    if (!value %in% c("free","series","grid","none"))
        stop("Invalid permutation type")
    value <- rep(value, length.out = 1L)
    object[["type"]] <- value
    object <- fixupCall(object, "type", value)
    object
}

## setMirror ############################################################
`setMirror<-` <- function(object, value) {
    UseMethod("setMirror<-")
}

`setMirror<-.default` <- function(object, value) {
    stop("No default method for `setMirror`")
}

`setMirror<-.how` <- function(object, value) {
    stop("`setMirror` can not be used directly on '\"how\"' objects.")
}

`setMirror<-.Within` <- function(object, value) {
    if (!is.null(value))
        value <- rep(as.logical(value), length.out = 1)
    object[["mirror"]] <- value
    object <- fixupCall(object, "mirror", value)
    object
}

`setMirror<-.Plots` <- function(object, value) {
    if (!is.null(value))
        value <- rep(as.logical(value), length.out = 1)
    object[["mirror"]] <- value
    object <- fixupCall(object, "mirror", value)
    object
}

## setConstant ############################################################
`setConstant<-` <- function(object, value) {
    UseMethod("setConstant<-")
}

`setConstant<-.default` <- function(object, value) {
    stop("No default method for `setConstant`")
}

`setConstant<-.how` <- function(object, value) {
    stop("`setConstant` can not be used directly on '\"how\"' objects.")
}

`setConstant<-.Within` <- function(object, value) {
    if (!is.null(value))
        value <- rep(as.logical(value), length.out = 1)
    object[["constant"]] <- value
    object <- fixupCall(object, "constant", value)
    object
}

`setConstant<-.Plots` <- function(object, value) {
    stop("`setConstant` does not apply to '\"Plots\"' objects.")
}
