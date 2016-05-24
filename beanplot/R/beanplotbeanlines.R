`beanplotbeanlines` <-
function (groups, side, beanlines, beanlinewd, at, boxwex, n, 
    col, horizontal, mlog, mexp) 
{
	res = NULL
    if (side == 4) 
        at <- rep(at, each = 2)
    for (i in 1:n) {
        if (beanlines == 2) {
            quants <- mexp(median(mlog(groups[[i]])))
        }
        else if (beanlines == 3) {
            quants <- mexp(quantile(mlog(groups[[i]])))
        }
        else {
            quants <- mexp(mean(mlog(groups[[i]])))
        }
		res = c(res, quants)
        x1 <- if (side == 3 || ((side == 4) && (i%%2 == 0))) 
            at[i]
        else at[i] - boxwex/2
        x2 <- if (side == 2 || ((side == 4) && (i%%2 == 1))) 
            at[i]
        else at[i] + boxwex/2
        if (horizontal) 
            segments(quants, x1, quants, x2, lwd = beanlinewd, 
                col = col[[i]][4])
        else segments(x1, quants, x2, quants, lwd = beanlinewd, 
            col = col[[i]][4])
    }
	res
}

