##' Plot a test function in 2D.
##'
##' @param x Function to plot.
##' @param lower Lower bounds of x1 and x2.
##' @param upper Upper bounds of x1 and x2.
##' @param n Number of locations at which to sample the function.
##' @param xlab Label of x (x1) axes.
##' @param ylab Label of y (x2) axes.
##' @param main Main title of plot.
##' @param log If \code{TRUE}, the z axes is plotted on log scale.
##' @param rank If \code{TRUE}, instead of the y values, their ranks are drawn.
##' @param show A vector of parts to plot. Defaults to
##'   \code{c("image", "contour")} and can be any subset.
##' @param ... Ignored.
##' @param image_args List of further arguments passed to image().
##' @param contour_args List of further arguments passed to contour().
##'
##' @examples
##' par(mfrow=c(2, 2))
##' f <- sphere_function(2)
##' plot(f)
##' plot(f, show="contour")
##' plot(f, rank=TRUE)
##' plot(f, log=TRUE)
##' @author Olaf Mersmann \email{olafm@@datensplitter.net}
##' @S3method plot soo_function
##' @method plot soo_function
##' @export
plot.soo_function <- function(x,
                              lower=lower_bounds(x),
                              upper=upper_bounds(x),
                              n=10000L,
                              main=function_name(x),
                              xlab=expression(x[1]),
                              ylab=expression(x[2]),
                              log=FALSE, rank=FALSE,
                              show=c("image", "contour"),
                              ...,
                              image_args=list(useRaster=TRUE),
                              contour_args=list()) {
  
  stopifnot(number_of_parameters(x) == 2,
            is.list(image_args),
            is.list(contour_args),
            n == as.integer(n))
  breaks_per_axis <- floor(sqrt(n))
  x1 <- seq(lower[1], upper[1], length.out=breaks_per_axis)
  x2 <- seq(lower[2], upper[2], length.out=breaks_per_axis)
  X <- expand.grid(x1, x2)
  z <- apply(X, 1, x)

  ## Shoud we logarithm the function values?
  if (log) {
    ## Ooops, some below zero. Lets fix that.
    if (any(z < 0)) {
      warning("Negative function values. Shifting function to apply logarithm.")
      z <- z - min(z) + 1
    }
    z <- log(z)
  }

  if (rank)
    z <- rank(z)

  ## Make z a breaks_per_axis times breaks_per_axis matrix:
  dim(z) <- c(breaks_per_axis, breaks_per_axis)

  ## Fixup image_args:
  if (!"col" %in% names(image_args))
    image_args$col <- terrain.colors(255)

  ## Plot image:
  if ("image" %in% show) {
    image_args <- append(list(x=x1, y=x2, z=z, xlab=xlab, ylab=ylab, main=main),
                         image_args)
    do.call(image, image_args)
  }
  ## Plot contour:
  if ("contour" %in% show) {    
    contour_args <- append(list(x=x1, y=x2, z=z), contour_args)
    ## If we did not plot an image, fixup contour args to include axis
    ## labels and plot title:
    if (!"image" %in% show)
      contour_args <- append(contour_args,
                             list(xlab=xlab, ylab=ylab, main=main))
    else
      contour_args$add <- TRUE
    do.call(contour, contour_args)
  }
}
