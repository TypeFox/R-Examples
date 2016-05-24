distfree.est <- function (n = NULL, alpha = NULL, P = NULL, side = 1) 
{
    dist.free.est2 <- function(n = NULL, alpha = NULL, P = NULL, 
        side = 1) {
        temp <- sum(c(is.null(n), is.null(alpha), is.null(P)))
        if (temp > 1) {
            stop(paste("Must specify values for any two of n, alpha, and P!", 
                "\n"))
        }
        if (side != 1 & side != 2) {
            stop(paste("Must specify a 1-sided or 2-sided interval!", 
                "\n"))
        }
        if (side == 1) {
            f <- function(n, alpha, P) alpha - P^n
            if (is.null(n)) {
                ret <- ceiling(uniroot(f, interval = c(0, 1e+100), 
                  alpha = alpha, P = P)$root)
            }
            if (is.null(P)) {
                ret <- ceiling(uniroot(f, interval = c(0, 1), 
                  n = n, alpha = alpha)$root * 10000)/10000
            }
            if (is.null(alpha)) {
                ret <- ceiling((1 - uniroot(f, interval = c(0, 
                  1), n = n, P = P)$root) * 10000)/10000
            }
        }
        else {
            f <- function(n, alpha, P) alpha - (n * P^(n - 1) - 
                (n - 1) * P^n)
            if (is.null(n)) {
                ret <- ceiling(uniroot(f, interval = c(0, 1e+100), 
                  alpha = alpha, P = P)$root)
            }
            if (is.null(P)) {
                ret <- ceiling(uniroot(f, interval = c(0, 1), 
                  n = n, alpha = alpha)$root * 10000)/10000
            }
            if (is.null(alpha)) {
                ret <- ceiling((1 - uniroot(f, interval = c(0, 
                  1), n = n, P = P)$root) * 10000)/10000
            }
        }
        ret
    }
    if (is.null(n)) {
        A <- length(alpha)
        B <- length(P)
        out <- matrix(0, nrow = A, ncol = B)
        for (i in 1:A) {
            for (j in 1:B) {
                out[i, j] <- dist.free.est2(alpha = alpha[i], 
                  P = P[j], side = side)
            }
        }
        rownames(out) <- alpha
        colnames(out) <- P
    }
    if (is.null(alpha)) {
        A <- length(n)
        B <- length(P)
        out <- matrix(0, nrow = A, ncol = B)
        for (i in 1:A) {
            for (j in 1:B) {
                out[i, j] <- dist.free.est2(n = n[i], P = P[j], 
                  side = side)
            }
        }
        rownames(out) <- n
        colnames(out) <- P
    }
    if (is.null(P)) {
        A <- length(alpha)
        B <- length(n)
        out <- matrix(0, nrow = A, ncol = B)
        for (i in 1:A) {
            for (j in 1:B) {
                out[i, j] <- dist.free.est2(alpha = alpha[i], 
                  n = n[j], side = side)
            }
        }
        rownames(out) <- alpha
        colnames(out) <- n
    }
    out
}
