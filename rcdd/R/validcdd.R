validcdd <- function(x, representation = c("H", "V")) {

    representation <- match.arg(representation)
    fred <- attr(x, "representation")
    if (! is.null(fred))
        representation <- match.arg(fred, c("H", "V"))

    if (! is.matrix(x))
        stop("ccd object must be matrix")

    if (! (is.character(x) | is.numeric(x)))
        stop("cdd object must be character or numeric")
    if (is.character(x))
        x <- q2q(x)

    if (ncol(x) <= 2)
        stop("cdd object must have more than two columns")

    if (! all(is.element(x[ , 1], 0:1)))
        stop("column one of ccd object must be zero-or-one valued")

    if (representation == "V")
        if (! all(is.element(x[ , 2], 0:1)))
            stop("column two of V-representation must be zero-or-one valued")

    return(TRUE)
}
