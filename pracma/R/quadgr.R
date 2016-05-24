##
##  q u a d g r . R  Gauss-Richardson Quadrature
##


quadgr <- function(f, a, b, tol = .Machine$double.eps^(1/2), ...) {
    stopifnot(is.numeric(a), length(a) == 1,
              is.numeric(b), length(b) == 1)

    fun <- match.fun(f)
    f <- function(x) fun(x, ...)

    # check order of limits
    if (a == b) return(list(value = 0, rel.err = 0))
    else if (a > b) {
        tmp <- a; a <- b; b <- tmp
        rev_p <- TRUE
    } else
        rev_p <- FALSE

    # Check infinite limits
    if (is.infinite(a) || is.infinite(b)) {
        if (is.finite(a) && is.infinite(b)) {
            f1 <- function(t) f(a + t/(1-t)) / (1-t)^2
            Q <- quadgr(f1, 0, 1, tol = tol)
    
        } else if (is.infinite(a) && is.finite(b)) {
            f2 <- function(t) f(b + t/(1+t)) / (1+t)^2
            Q <- quadgr(f2, -1, 0, tol = tol)
    
        } else if (is.infinite(a) && is.infinite(b)) {
            f1 <- function(t) f(t/(1-t)) / (1-t)^2
            f2 <- function(t) f(t/(1+t)) / (1+t)^2
            Q1 <- quadgr(f1,  0, 1, tol = tol/2)
            Q2 <- quadgr(f2, -1, 0, tol = tol/2)
            Q  <- list(value = Q1$value + Q2$value,
                       rel.err = Q1$rel.err + Q2$rel.err)
        }
        if (rev_p) Q$value <- -Q$value
        return(Q)
    } # else ...

    # 12-point Gauss-Legendre quadrature
    xq <- c(0.12523340851146894, 0.36783149899818018, 0.58731795428661748,
            0.76990267419430469, 0.9041172563704748,  0.98156063424671924)
    wq <- c(0.24914704581340288, 0.23349253653835478, 0.20316742672306584,
            0.16007832854334636, 0.10693932599531818, 0.047175336386511842)
    xq <- matrix(c(xq, -xq), ncol = 1)
    wq <- c(wq,  wq)
    nq <- length(xq)

    # Initiate vectors
    maxit <- 17                         # max number of iterations
    Q0 <- zeros(maxit,1)       	        # quadrature
    Q1 <- zeros(maxit,1)       	        # first Richardson extrapolation
    Q2 <- zeros(maxit,1)       	        # second Richardson extrapolation

    # One interval
    hh <- (b - a)/2                     # half interval length
    x <- (a + b)/2 + hh*xq              # nodes
    Q0[1] = hh * wq %*% fun(x)          # quadrature

    for (k in 3:maxit) {
        hh <- hh/2
        x <- cbind(x + a, x + b) / 2
        Q0[k] <- hh * wq %*% apply(f(x), 1, sum)

        # Richardson extrapolation
        if (k >= 5) {
            Q1[k] <- .rich(Q0,k)
            Q2[k] <- .rich(Q1,k)
        } else if (k >= 3) {
            Q1[k] <- .rich(Q0,k)
        }

        # Estimate absolute error
        if (k >= 6) {
            Qv <- c(Q0[k], Q1[k], Q2[k])
            Qw <- c(Q0[k-1], Q1[k-1], Q2[k-1])
        } else if (k >= 4) {
            Qv <- c(Q0[k], Q1[k])
            Qw <- c(Q0[k-1], Q1[k-1])
        } else {
            Qv <- Q0[k]
            Qw <- Q0[k-1]
        }

        err <- min(abs(Qv - Qw))
        j <- which.min(abs(Qv - Qw))
        Q <- Qv[j]
        
        # Convergence
        if (err < tol || !is.finite(Q)) break
    }

    if (rev_p) Q <- -Q
    return(list(value = Q, rel.err = err))
}


.rich <- function(Q, k) {
    if (Q[k] != Q[k-1]) {
        cc <- (Q[k-1] - Q[k-2]) / (Q[k] - Q[k-1]) - 1
    } else {
        cc <- 1
    }
    cc <- max(cc, 0.07)
    return(Q[k] + (Q[k] - Q[k-1])/cc)
}
