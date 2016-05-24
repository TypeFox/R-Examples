divchain.default <- function (x,y,z,...) {
#
    if(missing(z)) {
        if(inherits(x,"ppp") && is.factor(x$marks)) {
            z <- x$marks
        } else {
            stop("Either argument \"z\" was not supplied or it was not a factor.\n")
        }
    } else if(!is.factor(z)) stop("Argument \"z\" must be a factor.\n")
    dd <- deldir(x,y,z=z,...)
    divchain(dd)
}
