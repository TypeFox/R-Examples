
## this is definitely insane
##
## ceeboo 2007, 2008, 2014, 2015

setGeneric("format")

.list2S4 <- function(from, to) {
    if (!length(from))
        return (new(to))
    p <- pmatch(names(from), slotNames(to))
    if (any(is.na(p)))
        stop(paste("invalid slot name(s) for class", to, ":",
             paste(names(from)[is.na(p)], collapse = "")))
    names(from) <- slotNames(to)[p]
    do.call("new", c(from, Class = to))
}

.S42list <- function(object) {
    object <- attributes(object)
    object$class <- NULL
    object
}

setClass("SPparameter",
    representation(
        support = "numeric",
        maxsize = "integer",
        maxlen  = "integer",
        mingap  = "integer",
        maxgap  = "integer",
        maxwin  = "integer"
    ),

    prototype(support = 0.1, 
              maxsize = 10L, 
              maxlen  = 10L),

    validity = function(object) {
        if (object@support < 0 || object@support > 1)
            return("slot support : invalid range")
        if (length(object@maxsize) && object@maxsize < 1)
            return("slot maxsize : invalid range")
        if (length(object@maxlen) && object@maxlen < 1)
            return("slot maxlen : invalid range")
        if (length(object@mingap) && object@mingap < 1)
            return("slot mingap : invalid range")
        if (length(object@maxgap) && object@maxgap < 0)
            return("slot maxgap : invalid range")
        if (length(object@maxwin) && object@maxwin < 0)
            return("slot maxwin : invalid range")

        TRUE
    }
)

setMethod("initialize", "SPparameter",
    function(.Object, support, ...) {
        if (!missing(support))
            .Object@support <- support
        args <- list(...)
        for (name in names(args))
            slot(.Object, name) <- 
                as(args[[name]], Class = class(slot(.Object, name)))
        validObject(.Object)
        .Object
    }
)

setAs("NULL", "SPparameter",
    function(from, to) new(to))
setAs("list", "SPparameter", 
    function(from, to) .list2S4(from , to))
setAs("SPparameter", "list", 
    function(from) .S42list(from))
setAs("SPparameter", "vector",
    function(from) unlist(as(from, "list")))
setAs("SPparameter", "character", 
    function(from, to) unlist(lapply(as(from, "list"), as, class(to))))

setAs("SPparameter", "data.frame",
    function(from) {
        from <- as(from, "character")
        data.frame(name  = names(from), 
                   value = from, row.names = seq(from))
    }
)

.formatSP <- 
function(x, ...) {
    x <- as(x, "character")
    paste(format(names(x)), format(x, justify = "right"), sep = " : ")
}

setMethod("format", "SPparameter", 
    .formatSP
)

setMethod("show", signature(object = "SPparameter"),
    function(object) {
        out <- .formatSP(object)
        cat("set of", length(out), "spade parameters\n\n")
        if (length(out)) 
            cat(out, sep = "\n")

        invisible(NULL)
    }
)

## notes: (1) for now, we do not provide the BFS switch.
##        (2) the implicit default memory size is 32MB.

setClass("SPcontrol",
    representation(
        memsize = "integer",
        numpart = "integer",
        bfstype = "logical",
        verbose = "logical",
        summary = "logical",
       tidLists = "logical"
    ),

    prototype(bfstype = FALSE, verbose = FALSE, summary = FALSE,
	     tidLists = FALSE),

    validity = function(object) {
        if (length(object@memsize) && object@memsize < 16)
            return("slot memsize : invalid range")
        if (length(object@numpart) && object@numpart < 1)
            return("slot numpart : invalid range")

        TRUE
    }
)

setMethod("initialize", "SPcontrol",
    function(.Object, ...) {
        args <- list(...)
        for (name in names(args))
            slot(.Object, name) <-
                as(args[[name]], Class = class(slot(.Object, name)))
        validObject(.Object)
        .Object
    }
)

setAs("NULL", "SPcontrol",
    function(from, to) new(to))
setAs("list", "SPcontrol",
    function(from, to) .list2S4(from , to))
setAs("SPcontrol", "list",
    function(from) .S42list(from))
setAs("SPcontrol", "vector",
    function(from) unlist(as(from, "list")))
setAs("SPcontrol", "character", 
    function(from, to) unlist(lapply(as(from, "list"), as, class(to))))

setAs("SPcontrol", "data.frame",
    function(from) {
        from <- as(from, "character")
        data.frame(name  = names(from), 
                   value = from, row.names = seq(from))
    }
)

setMethod("format", "SPcontrol",
    .formatSP
)

setMethod("show", signature(object = "SPcontrol"),
    function(object) {
        out <- .formatSP(object)
        cat("set of", length(out), "spade control parameters\n\n")
        if (length(out))
            cat(out, sep = "\n")

        invisible(NULL)
    }
)

###
