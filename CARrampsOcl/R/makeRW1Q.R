makeRW1Q <-
function(n) {
        Q <- matrix(0, nrow = n, ncol = n)
        Q[1, 1:2] <- Q[n, n:(n-1)] <- c(1,-1)
        for( j in 2:(n-1) )
              Q[j,] <- c( rep(0,j-2), c(-1,2,-1), rep(0, n-j-1))
        Q
}

