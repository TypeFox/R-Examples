##' @title Extract Gradient Locations
##' @description Extract the gradient locations at which response curves were evaluated or for which counts were simulated.
##' @param x an object with \code{locations} as an attribute or a component, such as the object returned by \code{\link{coenocline}}.
##' @param ... arguments passed to other methods.
##' @return A vector or a matrix of gradient locations. For single-gradient simulations, a vector is returned, whereas for two-gradient simulations, a matrix of location pairs is returned.
##' @author Gavin L. Simpson
##'
##' @export
##'
##' @rdname locations
##'
##' @examples
##'
##' ## Poisson counts along a single gradient, Gaussian response
##' ## =========================================================
##'
##' x <- seq(from = 4, to = 6, length = 100)
##' opt <- c(3.75, 4, 4.55, 5, 5.5) + 0.5
##' tol <- rep(0.25, 5)
##' h <- rep(20, 5)
##'
##' ## simulate
##' set.seed(1)
##' y <- coenocline(x, responseModel = "gaussian",
##'                 params = cbind(opt = opt, tol = tol, h = h),
##'                 countModel = "poisson")
##' head(locations(y))
`locations` <- function(x, ...) {
    UseMethod("locations", x)
}

##' @rdname locations
##' @export
`locations.default` <- function(x, ...) {
    locs <- if (!is.null(locs <- attr(x, "locations"))) {
        locs
    } else if (!is.null(locs <- x[["locations"]])) {
        locs
    } else {
        stop("Couldn't find any gradient locations.")
    }
    locs
}
