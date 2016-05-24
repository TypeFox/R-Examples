linearity <- function(input, representation = c("H", "V")) {

    representation <- match.arg(representation)
    fred <- attr(input, "representation")
    if (! is.null(fred))
        representation <- match.arg(fred, c("H", "V"))
    h <- representation == "H"

    validcdd(input, representation)

    if (is.character(input)) {
        out <- .Call("impliedLinearity", input, h, PACKAGE = "rcdd")
    } else {
        if (! is.numeric(input))
            stop("input must be numeric or character")
        storage.mode(input) <- "double"
        out <- .Call("impliedLinearity_f", input, h, PACKAGE = "rcdd")
    }

    return(out)
}
