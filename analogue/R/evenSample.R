##' @title Number of samples per gradient segments
##'
##' @description The number of samples in sections along the gradient
##' is a useful diagnostic as to the quality of reconstructions at
##' gradient values within those sections.
##'
##' @details The sampling design of a training set, i.e. the number of
##' samples taken at points along the gradient, can influence the
##' uncertainty in the transfer function predictions at those values
##' of the gradient. Poorly sampled sections of the gradient may have
##' far larger RMSEP than the overall model RMSEP.
##'
##' @param grad numeric; vector of gradient values
##' @param n numeric; number of segments to partition the gradient into
##'
##' @return Numeric vector of length \code{n} containing the numbers of
##' samples per gradient segment.
##'
##' @author Gavin L. Simpson
##'
##' @keywords utilities
##'
##' @examples
##' data(SumSST)
##' ev <- evenSample(SumSST)  ## not an even sample...
##' ev
##' plot(ev)
`evenSample` <- function(grad, n = 10) {
    r <- range(grad)
    d <- diff(r)
    brks <- seq.int(r[1L], r[2L], length.out = n + 1)
    r <- r + (c(-1,1) * (d/1000))
    brks[c(1L, n + 1)] <- r
    segs <- findInterval(grad, brks)
    H <- hist(grad, breaks = brks, plot = FALSE)
    Nseg <- H$counts
    attr(Nseg, "gradient") <- deparse(substitute(grad))
    attr(Nseg, "numSegments") <- n
    attr(Nseg, "breaks") <- brks
    class(Nseg) <- "evenSample"
    Nseg
}


`print.evenSample` <- function(x, digits = getOption("digits"), ...) {
    attrs <- attributes(x)
    attributes(x) <- NULL
    x <- unclass(x)
    names(x) <- seq_along(x)
    writeLines(strwrap(paste("Gradient:", attrs[["gradient"]]), prefix = "\n"))
    writeLines(strwrap(paste("Segments:", attrs[["numSegments"]])))
    writeLines(strwrap("Number of samples per segment:", prefix = "\n"),
               sep = "\n\n")
    xx <- format(x, digits = digits, justify = "none")
    print(xx, quote = FALSE, right = TRUE, ...)
    invisible(x)
}
