getextent <-
function (x, y, ext) 
{
    stopifnot(is.numeric(x))
    stopifnot(is.numeric(y))
    if (class(ext)[1] == "Extent") {
        mnx <- xmin(ext)
        mxx <- xmax(ext)
        mny <- ymin(ext)
        mxy <- ymax(ext)
        xlm <- c(mnx, mxx)
        ylm <- c(mny, mxy)
    }
    if (is.raster(ext)) {
        mnx <- xmin(ext)
        mxx <- xmax(ext)
        mny <- ymin(ext)
        mxy <- ymax(ext)
        xlm <- c(mnx, mxx)
        ylm <- c(mny, mxy)
    }
    if (length(ext) == 4) {
        xlm <- ext[1:2]
        ylm <- ext[3:4]
    }
    if (length(ext) == 1) {
        mnx <- min(x)
        mxx <- max(x)
        mny <- min(y)
        mxy <- max(y)
        rgx <- (mxx - mnx) * 0.1
        rgy <- (mxy - mny) * 0.1
        xlm <- c(mnx - rgx, mxx + rgx)
        ylm <- c(mny - rgy, mxy + rgy)
    }
    fx1 <- (x < xlm[1]) * 1
    fx2 <- (x > xlm[2]) * 1
    fy1 <- (y < ylm[1]) * 1
    fy2 <- (y > ylm[2]) * 1
    beyond <- (fx1 + fx2 + fy1 + fy2 >= 1) * 1
    ex <- list(xlm = xlm, ylm = ylm, beyond = beyond)
    return(ex)
}
