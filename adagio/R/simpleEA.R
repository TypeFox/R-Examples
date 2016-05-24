##
##  s i m p l e E A . R  Simple Evolutionary Algorithm
##


simpleEA <-
function(fn, lower, upper, N = 100, ..., con = 0.1, new = 0.05,
         tol = 1e-10, eps = 1e-7, scl = 1/2, confined = FALSE, log = FALSE)
{
    stopifnot(is.numeric(lower), is.numeric(upper),
              length(lower) > 0, length(upper) > 0)
    if ( length(lower) != length(upper) ) 
        stop("Length of 'lower' and 'upper' limits must be the same.")
    if ( any(lower >= upper) )
        stop("Boundary: 'lower' < 'upper' not satisfied for all entries.")

    stopifnot(is.numeric(N), is.numeric(con), is.numeric(new),
              is.numeric(tol), is.numeric(eps), is.numeric(scl))
    if ( length(scl) != 1 || scl >= 1 || scl <= 0)
        stop("Skalar 'scl' must be positive and smaller than one.")

    fn  <- match.fun(fn)
    if (confined) {
        fun <- function(x) {
            if (any(x < lower) || any(x > upper)) Inf
            else fn(x, ...)
        }
    } else {
        fun <- function(x) fn(x, ...)
    }

    N <- floor(N)
    if (N <= 20)
        stop("No. of individuums per generation should be 20 or greater.")
    M <- floor(con * N)         # no. of individuals near/per parent
    K <- floor(new * N)         # no. of new/independent individuals

    if (M < 1 && K < N+1)
        stop("Argument 'con' too small: 'con > 1/N' must be true.")
    if (K < 0)      K <- 0
    if (M + K < N)  K <- N - M

    L <- c(Inf, Inf, Inf)

    n  <- length(lower)
    z0 <- (upper + lower)/2     # midpoint of rectangle
    h0 <- upper - lower         # side lengths of rectangle
    h <- scl * h0

    Parents <- rbind(z0, matrix(runif((N-1)*n), ncol=n) %*% diag(h0))
    fvals <- apply(Parents, 1, fun)

    niter <- 0
    while ( min(h) > eps ) {
        newParents <- rbind(
            Parents,
            kronecker(Parents, rep(1, M)) + 
                matrix(runif(M*N*n, -1, 1), ncol=n) %*% diag(h),
            matrix(runif(K*n), ncol=n) %*% diag(h) )

        # newFvals <- apply(newParents, 1, fun)
        newFvals <- c(fvals, apply(newParents[(N+1):(M*N+K), ], 1, fun))
        niter <- niter + 1

        new_order <- order(newFvals)[1:N]
        Parents <- newParents[new_order, ]
        fvals <- newFvals[new_order]

        newmin <- min(fvals)
        if (log) cat(niter, "\t", newmin, "\n")
        if (max(abs(L - newmin)) < tol) {
            break
        } else {
            L[2:3] <- L[1:2]
            L[1] <- newmin
        }

        # initialize new loop
        h <- scl * h
    }

    rel.tol <- max(abs(L - newmin))
    fun.calls <- N + niter * ((M-1) * N + K)
    return( list(par = Parents[1, ], val = min(fvals), fun.calls = fun.calls,
                 rel.scl = max(h), rel.tol = rel.tol) )
}


# simpleEA <-
# function(fun, lower, upper, N = 100, M = floor(0.1 * N),
#                      eps = 1e-6, scl = 1/2, log = FALSE)
# {
#     # if ( length(lower) != length(upper) ) error()
#     # if ( scale < 0 || scale >= 1 ) error
# 
#     n <- length(lower)
#     z0 <- (upper + lower)/2     # midpoint of rectangle
#     h0 <- upper - lower         # side lengths of rectangle
# 
#     Parents <- rbind(z0, matrix(runif((N-1)*n), ncol=n) %*% diag(h0))
#     fvals <- apply(Parents, 1, fun)
# 
#     h <- scl * h0
#     while ( min(h) > eps ) {
#         newParents <- rbind(
#             Parents,
#             kronecker(Parents, rep(1, M)) + 
#                 matrix(runif(M*N*n, -1, 1), ncol=n) %*% diag(h) )
# 
#         newFvals <- c(fvals, apply(newParents[(N+1):(M*N), ], 1, fun))
#         new_order <- order(newFvals)[1:N]
#         Parents <- newParents[new_order,]
#         fvals <- newFvals[new_order]
# 
#         if (log) cat(min(fvals), "\n")
# 
#         # initialize new loop
#         h <- scl * h
#     }
# 
#     return( list(par=Parents[1,], val=min(fvals), scl=max(h)) )
# }
