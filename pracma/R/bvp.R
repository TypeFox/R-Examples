##
##  b v p . R  Boundary Value Problems
##


bvp <- function(f, g, h, x, y, n = 50) {
    stopifnot(is.numeric(x), is.numeric(y), is.numeric(n))
    if (length(x) != 2 || length(y) != 2)
        stop("Arguments 'x' and 'y' must have length 2.")
    if (length(n) != 1 || floor(n) != ceiling(n) || n < 2)
        stop("Argument 'n' must be an integer greater or equal 2.")

    if (is.numeric(f)) ffun <- function(x) rep(f[1], length(x))
    else               ffun <- match.fun(f)
    if (is.numeric(g)) gfun <- function(x) rep(g[1], length(x))
    else               gfun <- match.fun(g)
    if (is.numeric(h)) hfun <- function(x) rep(h[1], length(x))
    else               hfun <- match.fun(h)

    xa <- x[1]; xb <- x[2]
    ya <- y[1]; yb <- y[2]
    xs <- linspace(xa, xb, n+2)[2:(n+1)]
    dt <- (xb - xa) / (n+1)

    a <- -2 - dt^2 * gfun(xs)           # main diagonal
    b <-  1 - dt/2 * ffun(xs[1:(n-1)])  # superdiagonal
    d <-  1 + dt/2 * ffun(xs[2:n])      # subdiagonal

    rhs <- dt^2 * hfun(xs)              # right hand side
    rhs[1] <- rhs[1] - ya * (1 + (dt/2) * ffun(xs[1]))
    rhs[n] <- rhs[n] - yb * (1 - (dt/2) * ffun(xs[n]))

    ys <- trisolve(a, b, d, rhs)
    return(list(xs = c(xa, xs, xb), ys = c(ya, ys, yb)))
}


# bvp <- function(f, a, b, ya, yb, N, cc, ...) {
#     stopifnot(is.numeric(a), length(a) == 1, is.numeric(b), length(b) == 1,
#               is.numeric(ya), length(ya) == 1, is.numeric(yb), length(yb) == 1)
# 
#     if (!is.numeric(N) || length(N) != 1 || floor(N) != ceiling(N) || N < 0)
#         stop("Argument 'N' must be an integer greater or equal 1.")
#     if (!is.numeric(cc) || length(cc) != 3)
#         stop("Argument 'cc' must be a real vector of length 3.")
# 
#     fun <- match.fun(f)
#     f   <- function(x) fun(x, ...)
# 
#     h <- (b-a)/(N+1)
#     xh <- linspace(a, b, N+2)
#     hm <- cc[1]/h^2
#     hd <- cc[2]/(2*h)
# 
#     A <- diag((2*hm + cc[3]), N, N)
#     A[col(A) == row(A)-1] <- -hm - hd
#     A[col(A) == row(A)+1] <- -hm - hd
# 
#     F <- f(xh[2:(N+1)])
#     F[1] <- F[1] + ya*(hm + hd)
#     F[N] <- F[N] + yb*(hm - hd)
# 
#     yh <- qr.solve(A, F)
#     yh <- c(ya, yh, yb)
# 
#     return(list(xh = xh, yh = yh))
# }
