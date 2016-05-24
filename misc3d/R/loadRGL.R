loadRGL <- function() {
    if (! suppressWarnings(require(rgl,quietly=TRUE)))
        stop("rgl is mot available")
}


