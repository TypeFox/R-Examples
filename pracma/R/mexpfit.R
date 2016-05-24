##
##  m e x p f i t . R  Multi-exponential Fitting
##


mexpfit <- function(x, y, p0, w = NULL, const = TRUE, options = list()) {
    stopifnot(is.numeric(x), is.numeric(y), is.numeric(p0))
    n <- length(x)
    if (length(y) != n)
        stop("Arguments 'x', 'y' must be of the same length.")
    p0 <- unique(p0)
    m <- length(p0)
    if (n <= 2*m+1)
        stop("Not enough data points available for fitting exponential sums.")

    opts <- list(tau     = 1e-4,
                 tolx    = 1e-7,
                 tolg    = 1e-9,
                 maxeval = 1000)
    namedOpts <- match.arg(names(options), choices = names(opts),
                           several.ok = TRUE)
    if (!is.null(names(options)))
        opts[namedOpts] <- options

    if (any(p0 == 0) || any(duplicated(p0)))
        stop("All entries in 'p0' must be different and not equal to zero.")

    .fexp <- function(b) {
        M <- outer(x, b, function(x, b) exp(b*x))
        if (const) M <- cbind(1, M)
        a <- qr.solve(M, y) 
        M %*% a - y
    }

    Lsq <- lsqnonlin(.fexp, p0, options = opts)
    b <- Lsq$x

    M <- outer(x, b, function(x, b) exp(b*x))
    if (const) M <- cbind(1, M)
    a <- qr.solve(M, y)

    if (const) { a0 <- a[1]; a <- a[-1] }
    else         a0 <- 0
    return(list(a0 = a0, a = a, b = b,
           ssq = Lsq$ssq, iter = Lsq$neval, errmess = Lsq$errmess))
}


lsqsep <- function(flist, p0, xdata, ydata, const = TRUE) {
    stopifnot(is.numeric(xdata), is.numeric(ydata), is.numeric(p0))
    n <- length(xdata)
    if (length(ydata) != n)
        stop("Numeric arguments 'xdata', 'ydata' must have the same length.")
    m <- length(flist)
    # lapply is.function

    .fapply <- function(b) {
        M <- matrix(1, nrow = n, ncol = m + 1)
        for (i in 1:m) {
            fi <- flist[[i]]
            xi <- fi(b[i], xdata)
            M[, i+1] <- xi
        }
        if (!const) M <- M[, 2:ncol(M)]
        a <- qr.solve(M, ydata)         # sum((M %*% a - ydata)^2)
        M %*% a - ydata                 # for lsqnonlin
    }

    # Find the function parameters b
    Lsq <- lsqnonlin(.fapply, p0)
    b <- Lsq$x

    # Find the linear parameters a
    M <- matrix(1, nrow = n, ncol = m + 1)
    for (i in 1:m) {
        fi <- flist[[i]]
        xi <- fi(b[i], xdata)
        M[, i+1] <- xi
    }
    if (!const) M <- M[, 2:ncol(M)]
    a <- qr.solve(M, ydata)

    if (const) {
        a0 <- a[1]; a <- a[2:length(a)]
    } else {
        a0 <- 0
    }
    return(list(a0 = a0, a = a, b = b, ssq = Lsq$ssq))
}
