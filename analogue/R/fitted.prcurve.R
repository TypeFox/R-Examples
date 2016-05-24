`fitted.prcurve` <- function(object, type = c("curve","smooths"), ...) {
    type <- match.arg(type)
    if (isTRUE(all.equal(type, "curve"))) {
        f <- object$s
    } else if (isTRUE(all.equal(type, "smooths"))) {
        f <- sapply(object$smooths, fitted)
        dimnames(f) <- dimnames(object$data)
        attr(f, "tag") <- object$tag
    } else {
        stop("Invalid 'type'")
    }
    f
}
