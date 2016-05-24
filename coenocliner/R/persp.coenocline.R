##' @title Perspective Plot of Species Simulations Along Gradients
##'
##' @description A simple S3 \code{\link{persp}} method for coenocline simulations.
##'
##' @param x an object of class \code{"coenocline"}, the result of a call to \code{\link{coenocline}}.
##' @param species vector indicating which species to plot. This can be any vector that you can use to subset a matrix, but numeric or logical vectors would be mostly commonly used.
##' @param theta,phi angles defining the viewing direction. \code{theta} gives the azimuthal direction and \code{phi} the colatitude. See \code{\link{persp}} for further details.
##' @param ... additional arguments to \code{\link{persp}}.
##'
##' @return A plot is drawn on the current device.
##'
##' @author Gavin L. Simpson
##'
##' @export
##'
##' @rdname persp.coenocline
##'
##' @keywords hplot
##'
##' @importFrom graphics persp
##'
##' @examples
##' ## Poisson counts along two correlated gradients, Gaussian response
##' ## ================================================================
##'
##' set.seed(1)
##' N <-  40
##' x1 <- seq(from = 4, to = 6, length = N)
##' opt1 <- seq(4, 6, length = 5)
##' tol1 <- rep(0.25, 5)
##' x2 <- seq(from = 2, to = 20, length = N)
##' opt2 <- seq(2, 20, length = 5)
##' tol2 <- rep(1, 5)
##' h <- rep(30, 5)
##' xy <- expand.grid(x = x1, y = x2)
##'
##' set.seed(1)
##' params <- list(px = list(opt = opt1, tol = tol1, h = h),
##'                py = list(opt = opt2, tol = tol2))
##' y <- coenocline(xy,
##'                 responseModel = "gaussian",
##'                 params = params,
##'                 extraParams = list(corr = 0.5),
##'                 countModel = "poisson")
##'
##' ## perspective plot(s) of simulated counts
##' layout(matrix(1:6, ncol = 3))
##' op <- par(mar = rep(1, 4))
##' persp(y)
##' par(op)
##' layout(1)
##'
##' ## as before but now just expectations
##' y <- coenocline(xy,
##'                 responseModel = "gaussian",
##'                 params = params,
##'                 extraParams = list(corr = 0.5),
##'                 countModel = "poisson",
##'                 expectation = TRUE)
##'
##' ## perspective plots of response curves
##' layout(matrix(1:6, ncol = 3))
##' op <- par(mar = rep(1, 4))
##' persp(y)
##' par(op)
##' layout(1)
##'
##' ## Same plots generated using the `plot` method
##' layout(matrix(1:6, ncol = 3))
##' op <- par(mar = rep(1, 4))
##' persp(y)
##' par(op)
##' layout(1)
`persp.coenocline` <- function(x, species = NULL, theta = 45, phi = 30, ...) {
    locs <- locations(x)
    if (NCOL(locs) != 2L) {
        stop("`persp()` method requires simulations along only two gradients.")
    } else {
        nspp <- NCOL(x)
        ## process species
        if (is.null(species) || missing(species)) {
            species <- seq_len(nspp)
        }
        if (length(species) > nspp)
            stop("'species' is longer than the number of species simulated")
    }

    ## order by first gradient, then second
    ord <- order(locs[,1], locs[,2])

    ## unique gradient locations
    locs1 <- unique(locs[,1])
    locs1 <- locs1[order(locs1)]
    locs2 <- unique(locs[,2])
    locs2 <- locs2[order(locs2)]
    nc <- length(locs2)

    ## loop to plot
    for (i in seq_along(species)) {
        sppi <- x[, i]                  # take the ith species to plot
        sppi <- matrix(sppi, ncol = nc) # reshape to a matrix
        persp(locs1, locs2, sppi, theta = theta, phi = phi, ...)  # plot using persp
    }
    invisible()
}
