##' perspEmpty
##'
##' Creates a persp plot without drawing anything in the interior.
##' Does equivalent \code{of plot( type="n")} for persp.
##'
##' Regression demonstrations require a blank slate in which
##' points and planes can be drawn. This function creates that
##' blank persp canvas for those projects. It is not necessary
##' that x1, x2 and y be vectors of the same length, since this
##' function's only purpose is to plot an empty box with ranges
##' determined by the input variables. persp calls the 3 axes
##' x, y, and z, but here they are called x1, x2, and y.
##'
##' @param x1 data for the first horizontal axis, an R vector
##' @param x2 data for the second horizontal axis, an R vector
##' @param y data for the vertical axis, an R vector
##' @param x1lab label for the x1 axis, (the one called "xlab" inside persp)
##' @param x2lab label for the x2 axis, (the one called "ylab" inside persp)
##' @param ylab label for the y (vertical) axis (the one called "zlab" inside persp)
##' @param x1lim Optional: limits for x1 axis (should be a vector with 2 elements)
##' @param x2lim Optional: limits for x2 axis (should be a vector with 2 elements)
##' @param ... further arguments that are passed to persp. Please note
##' Please remember that y is the vertical axis, but for persp, that
##' is the one I call x2.  Thus dot-dot-dot arguments including xlab,
##' ylab, zlab, xlim, ylim, and zlim are going to be ignored.
##' @return The perspective matrix that is returned by persp
##' @name perspEmpty
##' @export perspEmpty
##' @examples
##' x1 <- 1:10
##' x2 <- 41:50
##' y <-  rnorm(10)
##' perspEmpty(x1, x2, y)
##' res <- perspEmpty(x1, x2, y, ticktype="detailed", nticks=10)
##' mypoints1 <- trans3d ( x1, x2, y, pmat = res )
##' points( mypoints1, pch = 16, col= "blue")
perspEmpty <-
    function(x1, x2, y, x1lab = "x1", x2lab = "x2", ylab = "y", x1lim, x2lim, ... )
{
    x1range <- range(x1, na.rm = TRUE)
    x2range <- range(x2, na.rm = TRUE)
    yrange <- range(y, na.rm = TRUE)
    zZero <- outer(x1, x2, function(a,b) { a*b*0 + yrange[1] })

    dotargs <- list(...)
    dotargs[["xlab"]] <- x1lab
    dotargs[["ylab"]] <- x2lab
    dotargs[["zlab"]] <- ylab
    if (!missing(x1lim)) dotargs[["xlim"]] <- x1lim
    if (!missing(x2lim)) dotargs[["ylim"]] <- x2lim
    myDefaults <- list(x = x1, y = x2, z = zZero, zlim = yrange, lwd = 1, theta = -20, phi = 15)

    myargs <- modifyList(myDefaults, dotargs)
    res <- do.call("persp", myargs)
}

