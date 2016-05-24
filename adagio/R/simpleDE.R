##
##  s i m p l e D E . R  Simple Differential Evolution Algorithm
##


simpleDE <-
function(fun, lower, upper, N = 64, nmax = 256, r = 0.4, 
            confined = TRUE, log = FALSE)
{
    n <- length(lower)
    if (length(upper) != n)
        stop("'lower' and 'upper' must of the same length.")

    L <- matrix(rep(lower, each = N), nrow = N, ncol = n)
    U <- matrix(rep(upper, each = N), nrow = N, ncol = n)

    G <- matrix(runif(N*n), nrow = N, ncol = n)
    H <- G <- L + G * (U - L)
    F <- apply(G, 1, fun)

    for (g in 1:nmax) {
        for (i in 1:N) {
            ii <- sample(1:N, 3)
            ci <- G[ii[1], ] + r * (G[ii[2], ] - G[ii[3], ])

            if (confined) {
                if (any(ci < lower) || any(ci > upper)) {
                    ci <- lower + runif(n) * (upper - lower)
                }
            }

            fi <- fun(ci)
            if (fi < F[i]) {
                H[i, ] <- ci
                F[i] <- fi
            }
        }
        G <- H

        if (log && (g %% 10 == 0))
            cat(g, "  ", "\t", min(F), "\n", sep = "")
    }

    i0 <- which.min(F)
    return( list(fmin = F[i0], xmin = G[i0, ]) )
}
