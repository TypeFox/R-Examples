Mij <- function(x, s, i = 0, j = 0) {

    if(missing(s)) {

        xdim <- dim(x)
        s <- cbind(rep(1:xdim[1], xdim[2]), rep(1:xdim[2], each = xdim[1]))

    } # end of if missing 's' stmts.

    s[,1] <- s[,1]^i
    s[,2] <- s[,2]^j

    s <- apply(s, 1, prod, na.rm = TRUE)

    return(sum(s * c(x)))

} # end of 'Mij' function.

imomenter <- function(x, loc = NULL, ...) {

    UseMethod("imomenter", x)

} # end of 'imomenter' function.

imomenter.im <- function(x, loc = NULL, ...) {

    x <- as.matrix.im(x)
    UseMethod("imomenter", x)

} # end of 'imomenter.im' function.

imomenter.matrix <- function(x, loc = NULL, ...) {

    out <- list()
    # out$theCall <- match.call()

    if(is.null(loc)) {

        xdim <- dim(x)
        loc <- cbind(rep(1:xdim[1], xdim[2]), rep(1:xdim[2], each = xdim[1]))

    }

    M00 <- Mij(x, s = loc)
    M10 <- Mij(x, s = loc, i = 1)
    M01 <- Mij(x, s = loc, j = 1)
    M11 <- Mij(x, s = loc, i = 1, j = 1)
    M20 <- Mij(x, s = loc, i = 2)
    M02 <- Mij(x, s = loc, j = 2)

    xbar <- M10 / M00
    ybar <- M01 / M00

    cen <- c(xbar, ybar)
    names(cen) <- c("x", "y")

    mu11 <- M11 / M00 - xbar * ybar
    mu20 <- M20 / M00 - xbar^2
    mu02 <- M02 / M00 - ybar^2

    theta <- 0.5 * atan2(2 * mu11, mu20 - mu02)

    out$area <- M00
    out$centroid <- cen
    out$orientation.angle <- theta

    raw <- c(M00, M10, M01, M11, M20, M02)
    names(raw) <- c("M00", "M10", "M01", "M11", "M20", "M02")
    out$raw.moments <- raw

    out$cov <- rbind(c(mu20, mu11), c(mu11, mu02))

    class(out) <- "imomented"
    return(out)

} # end of 'imomenter.default' function.

print.imomented <- function(x, ...) {

    cat("Image Area: ", x$area, "\n")
    cat("Image centroid: \n")
    print(x$centroid)

    cat("Orientation Angle: ", x$orientation.angle, "\n")
    cat("Image covariance:\n")
    print(x$cov)

    cat("Raw moments:\n")
    print(x$raw.moments)

    invisible()

} # end of 'print.imomented' function.

