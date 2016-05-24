direct1d <- function(f, a, b, maxiter = 20, ...) {
    if (!is.numeric(a) || !is.numeric(b) ||
        length(a) != 1 || length(b) != 1)
        stop("Endpoints 'a', 'b' must be numeric and of length 1.")
    if (a >= b || !is.finite(a) || !is.finite(b))
        stop("Endpoints 'a', 'b' must be finite and a < b must hold.")

    fun <- match.fun(f)
    f <- function(x) fun(x, ...)

    #(1) initialize
    c0 <- (b-a)/2
    f0 <- f(c0); fcnt <- 1
    m <- 1
    B <- matrix(c(m, a, b, c0, f0), nrow = 1)
    xmin <- c0; fmin <- f0; tmin <- (b-a)/2
    niter <- 0

    while(niter < maxiter) {
        #(2) identify potentially optimal intervals
        s <- .identify_optimal(B)
        S <- B[s, , drop = FALSE]


        #(3) select all intervals j in S
        for (i in 1:nrow(S)) {
            k <- S[i, 1]

            #(4) divide into three parts
            d  <- (S[i, 3] - S[i, 2])/3
            c1 <- S[i, 4] - d; c2 <- S[i, 4] + d
            f1 <- f(c1); f2 <- f(c2); fcnt <- fcnt + 2
            # fmin <- min(f1, fmin, f2)
            if (f1 < fmin) {
                xmin <- c1; fmin <- f1; tmin <- d
            } else if (f2 < fmin) {
                xmin <- c2; fmin <- f2; tmin <- d
            }

            #(5) Add and update intervals
            m <- m+1; B <- rbind(B, c(m, S[i, 2],       S[i, 2] + d, c1, f1))
            m <- m+1; B <- rbind(B, c(m, S[i, 2] + 2*d, S[i, 3],     c2, f2))
            B[k, 2:3] <- c(S[i, 2] + d, S[i, 3] - d)
        }
        #(6) S is empty now
        
        #(7) check iteration limit
        niter <- niter + 1
     }
     return(list(xmin = xmin, fmin = fmin))
}


.identify_optimal <- function(B) {
    P <- cbind(B[, 1], (B[, 3] - B[, 2])/2, B[, 5])
    i0 <- 1
    for (i in 1:nrow(P)) {
        if ( P[i, 2] <  P[i0, 2]  ||
            (P[i, 2] == P[i0, 2]  && P[i, 3] < P[i0, 3]) ) i0 <- i
    }
    I <- P[i0, 1]
    p2 <- P[i0, 2]; p3 <- P[i0, 3]
    P <- P[P[, 2] > p2, , drop = FALSE]

    while (!isempty(P)) {
        t <- atan2(P[, 3, drop = FALSE] - p3, P[, 2] - p2)
        ti <- which.min(t)
        i0 <- P[ti, 1]
        I <- c(I, i0)
        p2 <- P[ti, 2]; p3 <- P[ti, 3]
        P <- P[P[, 2] > p2, , drop = FALSE]
    }
    return(I)
}
