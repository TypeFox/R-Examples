`beanplotinnerborders` <-
function (innerborder, at, dens, side, displayn, n, horizontal, 
    mexp) 
{
    if (side == 4) {
        x0 <- numeric(displayn)
        x1 <- numeric(displayn)
        for (i in 1:displayn) {
            if (i + 1 > n) {
                x0[i] <- dens[["x", i * 2 - 1]][1]
                x1[i] <- tail(dens[["x", i * 2 - 1]], 1)
            }
            else {
                x0[i] <- min(dens[["x", i * 2 - 1]][1], dens[["x", 
                  i * 2]][1])
                x1[i] <- max(tail(dens[["x", i * 2 - 1]], 1), 
                  tail(dens[["x", i * 2]], 1))
            }
        }
    }
    else {
        x0 <- unlist(lapply(dens["x", ], head, 1))
        x1 <- unlist(lapply(dens["x", ], tail, 1))
    }
    if (horizontal) 
        segments(mexp(x0), at, mexp(x1), at, col = innerborder)
    else segments(at, mexp(x0), at, mexp(x1), col = innerborder)
}

