### Extractor, replacement functions for blocks

### Extractor methods
`blocks` <- function(object, ...) {
    UseMethod("blocks")
}

`blocks.default` <- function(object, ...) {
    stop("No default method for 'blocks()'")
}

`blocks.how` <- function(object, ...) {
    object$blocks
}

### Replacement methods
`blocks<-` <- function(object, value) {
    UseMethod("setBlocks<-")
}

`blocks<-.default` <- function(object, value) {
    stop("No default method for `setBlocks`")
}

`blocks<-.how` <- function(object, value) {
    object[["blocks.name"]] <- deparse(substitute(value))
    if (!is.null(value))
        value <- as.factor(value)
    object["blocks"] <- list(value)
    object <- fixupCall(object, "blocks", value)
    object
}

### Temporarily drop this here
`Blocks` <- function(strata = NULL) {
    out <- list(strata = strata)
    ## keep as list for now
    ##class(out) <- "Blocks"
    out
}
