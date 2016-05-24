`beanplotscatters` <-
function (groups, side, method, jitter, dens, at, wd2, boxwex, 
    n, ll, maxstripline, col, horizontal, mlog, mexp) 
{
    if (side == 4) 
        at <- rep(at, each = 2)
    for (i in 1:n) {
        y <- groups[[i]]
        if (length(y) == 0) 
            next
        if (method == 2) {
            y <- aggregate(as.numeric(y), list(as.numeric(y)), 
                length)
            y[1] <- as.numeric(as.character(y[[1]]))
        }
        else if (method == 3) {
            y <- list(mexp(jitter(mlog(y), amount = jitter)), 
                rep(1, length(y)))
        }
        else {
            y <- list(y, rep(1, length(y)))
        }
        y[[2]] <- sapply(y[[2]] * ll * boxwex/2, min, maxstripline * 
            boxwex/2)
        if (col[[i]][2] != col[[i]][3]) {
            crossings <- approx(dens[["x", i]], dens[["y", i]], 
                mlog(y[[1]]))$y * wd2
            crossings <- pmin(crossings, y[[2]])
            extra <- unlist(y[2]) - crossings
        }
        else {
            crossings <- y[[2]]
            extra <- 0
        }
        if (max(extra) > 0) {
            w <- extra > 0
            if (horizontal) {
                if (side %in% c(1, 2) || ((side == 4) && (i%%2 == 
                  1))) 
                  segments(y[[1]][w], -extra[w] - crossings[w] + 
                    at[i], y[[1]][w], -crossings[w] + at[i], 
                    col = col[[i]][3])
                if (side %in% c(1, 3) || ((side == 4) && (i%%2 == 
                  0))) 
                  segments(y[[1]][w], crossings[w] + at[i], y[[1]][w], 
                    extra[w] + crossings[w] + at[i], col = col[[i]][3])
            }
            else {
                if (side %in% c(1, 2) || ((side == 4) && (i%%2 == 
                  1))) 
                  segments(-extra[w] - crossings[w] + at[i], 
                    y[[1]][w], -crossings[w] + at[i], y[[1]][w], 
                    col = col[[i]][3])
                if (side %in% c(1, 3) || ((side == 4) && (i%%2 == 
                  0))) 
                  segments(crossings[w] + at[i], y[[1]][w], extra[w] + 
                    crossings[w] + at[i], y[[1]][w], col = col[[i]][3])
            }
        }
        x1 <- if (side == 3 || ((side == 4) && (i%%2 == 0))) 
            at[i]
        else -crossings + at[i]
        x2 <- if (side == 2 || ((side == 4) && (i%%2 == 1))) 
            at[i]
        else crossings + at[i]
        if (horizontal) 
            segments(y[[1]], x1, y[[1]], x2, col = col[[i]][2])
        else segments(x1, y[[1]], x2, y[[1]], col = col[[i]][2])
    }
}

