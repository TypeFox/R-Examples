# draw a 3D cross (+) symbol in an rgl view


cross3d <- 
function(centre=rep(0,3), scale=rep(1,3), ...) {
    axes <- matrix(
      c(-1, 0, 0,   1, 0, 0,
        0, -1, 0,   0, 1, 0,
        0, 0, -1,   0, 0, 1),  6, 3, byrow=TRUE)
    if (!missing(scale)) {
        if (length(scale) != 3) scale <- rep(scale, length.out=3) 
        axes <- rgl::scale3d(axes, scale[1], scale[2], scale[3])
        }
    if (!missing(centre)) { 
        if (length(centre) != 3) scale <- rep(centre, length.out=3) 
        axes <- rgl::translate3d(axes, centre[1], centre[2], centre[3])
        }
    rgl::segments3d(axes, ...)
    invisible(axes)

}
