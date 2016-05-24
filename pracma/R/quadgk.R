##
##  g a u s s g k . R  Adapitve Gauss-Kronrod
##


quadgk <- function(f, a, b, tol = .Machine$double.eps^0.5, ...) {
    stopifnot(is.numeric(a), length(a) == 1,
              is.numeric(b), length(b) == 1)
    eps <- .Machine$double.eps

    fun <- match.fun(f)
    f <- function(x) fun(x, ...)
    
    if (a == b)     return(0)
    else if (a > b) return(-1 * quadgk(f, b, a, tol = tol))

    # Nodes and weights for Gauss-Kronrod (7, 15)
    n15 <- c(-0.9914553711208126, -0.9491079123427585, -0.8648644233597691,
             -0.7415311855993944, -0.5860872354676911, -0.4058451513773972,
             -0.2077849550078985,  0.0,                 0.2077849550078985,
              0.4058451513773972,  0.5860872354676911,  0.7415311855993944,
              0.8648644233597691,  0.9491079123427585,  0.9914553711208126)
    n7  <- c(-0.9491079123427585, -0.7415311855993944, -0.4058451513773972,
               0.0,
               0.4058451513773972, 0.7415311855993944,  0.9491079123427585)
    
    w15 <- c(0.02293532201052922, 0.06309209262997855,  0.1047900103222502,
             0.1406532597155259,  0.1690047266392679,   0.1903505780647854,
             0.2044329400752989,  0.2094821410847278,   0.2044329400752989,
             0.1903505780647854,  0.1690047266392679,   0.1406532597155259,
             0.1047900103222502,  0.06309209262997855,  0.02293532201052922)
    w7  <- c(0.1294849661688697,  0.2797053914892767,   0.3818300505051189,
             0.4179591836734694,
             0.3818300505051189,  0.2797053914892767,   0.1294849661688697)

    .gkadpt <- function(f, a, b, tol = tol) {
        # use nodes and weights from the environment
        x15 <- 0.5 * ((b - a) * n15 + b + a)
        x7  <- 0.5 * ((b - a) * n7  + b + a)
        Q7  <- sum(w7  * f(x7))  * (b-a)/2
        Q15 <- sum(w15 * f(x15)) * (b-a)/2

        if (!is.finite(Q7) || !is.finite(Q15)) {
            warning("Infinite or NA function value encountered.")
            return(Q15)
        } else if (abs(Q15 - Q7) < tol) {
            return(Q15)
        } else if (abs(b-a) < 16*eps) {
            warning("Minimum step size reached; singularity possible.")
            return(Q2)
        } # else

        Q2 <- .gkadpt(f, (a+b)/2, b, tol = tol)
        Q1 <- .gkadpt(f, a, (a+b)/2, tol = tol)

        return(Q1 + Q2)
    }

    # start the recursive procedure
    .gkadpt(f, a, b, tol = tol)
}
