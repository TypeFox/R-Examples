##
##  m a x s u b . R  Maximal Sum Subarray Problem
##  m a x s u b 2 d  Maximal Sum Subrectangle Problem
##


maxsub <- function(x, inds = TRUE, compiled = TRUE) {
    if (!is.numeric(x))
        stop("Argument 'x' must be a numeric vector.")
    n <- length(x)
    i1 <- 0; i2 <- 0
    s <- 0.0

    if (compiled) {
        S <- .Fortran("maxsubf", x = as.numeric(x), n = as.integer(n),
                                 s = as.numeric(s),
                                 i1 = as.integer(i1), i2 = as.integer(i2),
                                 PACKAGE = "adagio")
        if (inds)
            return(list(sum = S$s, inds = c(S$i1, S$i2)))
        else
            return(S$s)

    } else {
        if (!inds) {
            m1 <- m2 <- 0.0
            for (i in 1:n) {
                m2 <- max(m2 + x[i], 0.0)
                m1 <- max(m1, m2)
            }
            return(m1)
        } else {
            m1 <- m2 <- 0
            p1 <- p2 <- 0
            q1 <- q2 <- 1

            for (i in 1:n) {
                if (m2 > -x[i]) {
                    m2 <- m2 + x[i]
                    q2 <- i
                    if (m2 > m1) {
                        m1 <- m2
                        p1 <- q1; p2 <- q2
                    }
                } else {
                    m2 <- 0
                    q1 <- q2 <- i+1
                }
            }
            return(list(sum = m1, inds = c(p1, p2)))
        }
    }
}


maxsub2d <- function(A) {
    stopifnot(is.numeric(A), is.matrix(A))
    n <- nrow(A); m <- ncol(A)

    if (all(A <= 0))
        stop("At least on element of matrix 'A' must be positive.")
    if (all(A >= 0))
        return(list(sum = sum(A), inds = c(1, n, 1, m), submat = A))

    mi <- vector("integer", 4)
    S <- matrix(0, nrow = n+1, ncol = m)
    aa <- numeric(m)
    b <- 0.0

    fm <- 0.0
    R <- .Fortran("maxsub2f", as.numeric(A), as.numeric(S),
                              as.integer(n), as.integer(m),
                              fmax = as.numeric(fm), mind = as.integer(mi),
                              as.numeric(aa), as.numeric(b),
                              PACKAGE = "adagio")

    fm <- R$fmax
    mi <- R$mind

    return(list(sum = fm, inds = mi,
                submat = A[mi[1]:mi[2], mi[3]:mi[4], drop = FALSE]))
}


# maxsub2d <- function(A) {
#     stopifnot(is.numeric(A), is.matrix(A))
#     n <- nrow(A); m <- ncol(A)
# 
#     if (all(A <= 0))
#         stop("At least on element of matrix 'A' must be positive.")
#     if (all(A >= 0))
#         return(list(sum = sum(A), inds = c(1, n, 1, m), submat = A))
# 
#     S <- matrix(0, nrow = n+1, ncol = m)
#     S[1, ] <- 0
#     for (i in 2:(n+1))
#         S[i, ] <- S[i-1, ] + cumsum(A[i-1, ])
# 
#     mm <- 0
#     mi <- c(0, 0, 0, 0)
# 
#     for (i in 2:(n+1)) {
#         for (j in i:(n+1)) {
#             a <- numeric(n)
#             a[1] <- S[j, 1] - S[i-1, 1]
#             for (k in 2:n)
#                 a[k] <- S[j, k] - S[(i-1), k] - sum(a[1:(k-1)])
#             ms <- maxsub(a)
#             if (ms$sum > mm) {
#                 mm <- ms$sum
#                 mi <- c(i-1, j-1, ms$inds[1], ms$inds[2])
#             }
#         }
#     }
#     return(list(sum = mm, inds = mi,
#             submat = A[mi[1]:mi[2], mi[3]:mi[4], drop = FALSE]))
# }
