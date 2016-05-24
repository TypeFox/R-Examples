makeRW2Q <-
function(n) {
        Q <- matrix(0, nrow = n, ncol = n)
        Q[1, 1:3] <- Q[n, n:(n - 2)] <- c(1, -2, 1)
        Q[2, 1:4] <- Q[n - 1, n:(n - 3)] <- c(-2, 5, -4, 1)
        for (j in 3:(n - 2)) Q[j, ] <- c(rep(0, j - 3), c(1,
            -4, 6, -4, 1), rep(0, n - j - 2))
        Q
    }

