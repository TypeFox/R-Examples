`beanplotpolyshapes` <-
function (side, dens, at, wd2, combinedpolygons, displayn, n, 
    col, border, horizontal, mlog, mexp) 
{
    if ((side == 4) && (!combinedpolygons)) 
        at <- rep(at, each = 2)
    pborder <- NA
    for (i in 1:ifelse(combinedpolygons, displayn, n)) {
        if (combinedpolygons) {
            x1 <- rev(dens[["y", i * 2]] * wd2) + at[i]
            x2 <- dens[["y", i * 2 - 1]] * -wd2 + at[i]
            y1 <- rev(dens[["x", i * 2]])
            y2 <- dens[["x", i * 2 - 1]]
        }
        else {
            if (side == 2 || ((side == 4) && (i%%2 == 1))) {
                x1 <- NULL
                y1 <- NULL
            }
            else {
                x1 <- dens[["y", i]] * wd2 + at[i]
                y1 <- dens[["x", i]]
            }
            if (side == 3 || ((side == 4) && (i%%2 == 0))) {
                x2 <- NULL
                y2 <- NULL
            }
            else {
                x2 <- rev(dens[["y", i]]) * -wd2 + at[i]
                y2 <- rev(dens[["x", i]])
            }
        }
        if (length(x1) > 0) {
            x1 <- c(at[i], x1, at[i])
            y1 <- c(y1[1], y1, y1[length(y1)])
        }
        if (length(x2) > 0) {
            x2 <- c(at[i], x2, at[i])
            y2 <- c(y2[1], y2, y2[length(y2)])
        }
        if (combinedpolygons || (side == 1)) 
            pborder <- border[[i]]
        if (horizontal) {
            polygon(mexp(c(y1, y2)), c(x1, x2), col = col[[i]][1], 
                border = pborder)
            if ((!combinedpolygons) && (side != 1)) {
                if (is.null(border[[i]])) 
                  lines(mexp(c(y1, y2)), c(x1, x2))
                else if (!is.na(border[[i]])) 
                  lines(mexp(c(y1, y2)), c(x1, x2), col = border[[i]])
            }
        }
        else {
            polygon(c(x1, x2), mexp(c(y1, y2)), col = col[[i]][1], 
                border = pborder)
            if ((!combinedpolygons) && (side != 1)) {
                if (is.null(border[[i]])) 
                  lines(c(x1, x2), mexp(c(y1, y2)))
                else if (!is.na(border[[i]])) 
                  lines(c(x1, x2), mexp(c(y1, y2)), col = border[[i]])
            }
        }
    }
}

