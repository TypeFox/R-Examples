##
##  e l l i p . R  Elliptic Integrals
##


ellipke <- function(m, tol = .Machine$double.eps) {
    stopifnot(is.numeric(m))
    m <- c(m)
    if (any(m < 0) || any(m > 1))
        stop("Some elements of argument 'm' are out of range.")

    a0 <- 1
    b0 <- sqrt(1-m)
    s0 <- m
    i1 <- 0
    mm <- 1
    while (mm > tol) {
        a1 <- (a0+b0)/2
        b1 <- sqrt(a0*b0)
        c1 <- (a0-b0)/2
        i1 <- i1 + 1
        w1 <- 2^i1 * c1^2
        mm <- max(w1)
        s0 <- s0 + w1
        a0 <- a1
        b0 <- b1
    }
    k <- pi / (2*a1)
    e <- k * (1-s0/2)
    im <- finds(m == 1)
    if (!isempty(im)) {
        e[im] <- ones(length(im), 1)
        k[im] <- Inf
    }
    return(list(k = k, e = e))
}


##  Jacobi elliptic functions
ellipj <- function(u, m, tol = .Machine$double.eps) {
    stopifnot(is.numeric(u), is.numeric(m) || is.complex(m))
    u <- c(u); m <- c(m)
    if (length(u) == 1) {
        u <- rep(u, length(m))
    } else if (length(m) == 1) {
        m <- rep(m, length(u))
    } else {
        if (length(u) != length(m))
            stop("Arguments 'u' and 'm' must be of the same length.")
    }
    if (any(m < 0) || any(m > 1))
        stop("Some elements of argument 'm' are out of range.")

    mmax <- length(u); chunk <- 10
    cn <- sn <- dn <- numeric(mmax)
    a  <- b  <- cc <- matrix(0, nrow = chunk, ncol = mmax)
    a[1, ]  <- 1
    b[1, ]  <- sqrt(1-m)
    cc[1, ] <- sqrt(m)

    n <- numeric(mmax)
    i <- 1
    while (any(abs(cc[i, ]) > tol)) {
        i <- i + 1
        if (i > nrow(a)) {
          a  <- rbind(a,  matrix(0, chunk, mmax))
          b  <- rbind(b,  matrix(0, chunk, mmax))
          cc <- rbind(cc, matrix(0, chunk, mmax))
        }
        a[i, ]  <- 0.5 * (a[i-1, ] + b[i-1, ])
        b[i, ]  <-   sqrt(a[i-1, ] * b[i-1, ])
        cc[i, ] <- 0.5 * (a[i-1, ] - b[i-1, ])
        inds <- which(abs(cc[i, ]) <= tol & abs(cc[i-1, ]) > tol)
        if (!isempty(inds)) {
          mi <- 1; ni <- length(inds)  # [mi,ni] <- size(inds)
          n[inds]     <- rep(i-1, ni)  # repmat((i-1), mi, ni)
        }
    }

    phin <- matrix(0, nrow = i, ncol = mmax)
    phin[i, ] <- 2^n * a[i, ] * u
    while (i > 1) {
        i <- i - 1
        inds <- which(n >= i)
        phin[i, ] <- phin[i+1, ]
        if (!isempty(inds)) {
          phin[i, inds] <- 0.5 * (asin(cc[i+1, inds] * 
              sin(phin[i+1, inds] %% (2*pi)) / a[i+1, inds]) + phin[i+1,inds])
        }
    }
    # the general case
    sn <- sin(phin[1, ] %% (2*pi))
    cn <- cos(phin[1, ] %% (2*pi))
    dn <- sqrt(1 - m * sn^2)

    # some special cases
    m1 <- which(m == 1)         # special case m = 1
    sn[m1] <- tanh(u[m1])
    cn[m1] <- sech(u[m1])
    dn[m1] <- sech(u[m1])
    dn[m == 0] <- 1             # special case m = 0

    return(list(sn = sn, cn = cn, dn = dn))
}
