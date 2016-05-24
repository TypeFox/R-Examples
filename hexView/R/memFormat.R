
# A memFormat is a list of memBlock's  

is.memFormat <- function(x) {
    inherits(x, "memFormat")
}

memFormat <- function(...) {
    blocks <- list(...)
    if (!all(sapply(blocks, is.memBlock)))
        stop("Invalid format")

    class(blocks) <- "memFormat"
    blocks
}

