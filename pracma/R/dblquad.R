##
##  d b l q u a d . R  Double Integration
##


dblquad <- function(f, xa, xb, ya, yb, dim = 2, ..., 
                    subdivs = 300, tol = .Machine$double.eps^0.5) {
    stopifnot(is.numeric(xa), length(xa) == 1, is.numeric(xb), length(xb) == 1,
              is.numeric(ya), length(ya) == 1, is.numeric(yb), length(yb) == 1)

    fun <- match.fun(f)
    f <- function(x, y) fun(x, y, ...)
    if (length(f(c(xa, xb), c(ya, yb))) != 2)
        stop("Function 'f' does not appear to be vectorized.")

    if (dim == 2) {
        fy <- function(x)
                integrate(function(y) f(x, y), ya, yb,
                          subdivisions = subdivs, rel.tol = tol)$value
        Fy <- Vectorize(fy)
        Q  <- integrate(Fy, xa, xb,
                        subdivisions = subdivs, rel.tol = tol)$value

    } else if (dim == 1) {
        fx <- function(y)
                integrate(function(x) f(x, y), xa, xb,
                          subdivisions = subdivs, rel.tol = tol)$value
        Fx <- Vectorize(fx)
        Q  <- integrate(Fx, ya, yb,
                        subdivisions = subdivs, rel.tol = 10*tol)$value

    } else
        stop("Argument 'dim' can only be 1 (x-) or 2 (y-variable first).")

    return(Q)
}


triplequad <- function(f, xa, xb, ya, yb, za, zb, 
                        subdivs = 300, tol = .Machine$double.eps^0.5, ...) {
    stopifnot(is.numeric(xa), length(xa) == 1, is.numeric(xb), length(xb) == 1,
              is.numeric(ya), length(ya) == 1, is.numeric(yb), length(yb) == 1,
              is.numeric(za), length(za) == 1, is.numeric(zb), length(zb) == 1)

    fun <- match.fun(f)
    f <- function(x, y, z) fun(x, y, z, ...)

    fyz <- function(y, z) {
        Qin <- numeric(length(y))
        for (i in 1:length(y)) {
            fx  <- function(x) f(x, y[i], z[i])
            Qin <- integrate(fx, xa, xb,
                             subdivisions = subdivs, rel.tol = 1e-10)$value
        }
        Qin
    }
    fyz <- Vectorize(fyz)
    dblquad(fyz, ya, yb, za, zb, tol = tol)
}


simpson2d <- function(f, xa, xb, ya, yb, nx = 128, ny = 128, ...) {
    stopifnot(is.numeric(xa), length(xa) == 1, is.numeric(xb), length(xb) == 1,
              is.numeric(ya), length(ya) == 1, is.numeric(yb), length(yb) == 1)

    fun <- match.fun(f)
    f <- function(x, y) fun(x, y, ...)

    if (nx %% 2 != 0) nx <- nx + 1
    if (ny %% 2 != 0) ny <- ny + 1

    # Grid and grid vectors
    hx <- (xb - xa) / nx
    hy <- (yb - ya) / ny
    xg <- seq(xa, xb, by = hx)
    yg <- seq(ya, yb, by = hy)

    # Interchange meshgrid
    mgrid <- meshgrid(yg, xg)
    X <- mgrid$Y
    Y <- mgrid$X

    F <- f(X, Y)

    # Contributions from the corner points
    s1 <- F[1, 1] + F[1, ny+1] + F[nx+1, 1] + F[nx+1, ny+1]

    # Contributions from other edge points
    ixo <- seq(2, nx, by = 2); ixe <- seq(3, nx-1, by = 2)
    iyo <- seq(2, ny, by = 2); iye <- seq(3, ny-1, by = 2)
    s2 <- 2 * ( sum(F[1, iye]) + sum(F[nx+1, iye]) + sum(F[ixe, 1]) + sum(F[ixe, ny+1]) );
    s3 <- 4 * ( sum(F[1, iyo]) + sum(F[nx+1, iyo]) + sum(F[ixo, 1]) + sum(F[ixo, ny+1]) );

    # Contributions from interior points
    s4 <- 16 * sum( sum( F[ixo,iyo] ) ) + 4 * sum( sum( F[ixe,iye] ) );
    s5 <-  8 * sum( sum( F[ixe,iyo] ) ) + 8 * sum( sum( F[ixo,iye] ) );

    S <- hx * hy * (s1 + s2 + s3 + s4 + s5) / 9.0
    return(S)
}
