"dots" <-
function (x, y = 0.1, xlim = range(x, na.rm = TRUE), stacked = FALSE, 
    hmax = 0.5, base = TRUE, axes = FALSE, pch = 21, pch.size = "x", 
    labels = NULL, hcex = 1, cex = par("cex"), cex.axis = par("cex.axis")) 
{
    x <- x[!is.na(x)]
    hdots <- y
    xmin <- xlim[1]
    xmax <- xlim[2]
    x <- x[(x >= xmin) & (x <= xmax)]
    b <- strwidth(pch.size, units = "user", cex = cex)
    h <- strheight(pch.size, units = "user", cex = hcex * cex)
    if (stacked) {
        if (xmax - xmin < b) {
            stop("x-dimension resolution problem")
        }
        else {
            xu <- seq(xmin, xmax, by = b)
        }
        m <- length(xu)
        tab <- data.frame(j = 1:m, k = rep(0, m), xu = xu)
        n <- length(x)
        y <- rep(0, n)
        for (i in 1:n) {
            l <- max(tab$j[tab$xu <= x[i]])
            x[i] <- xu[l] + b/2
            tab$k[l] <- 1 + tab$k[l]
            y[i] <- tab$k[l]
        }
        y <- y * h
        u <- hdots + max(y)
        if (hmax <= hdots) 
            warning(paste("dot base <hdots=", hdots, "> higher than maximum column height <hmax=", 
                hmax, ">...", sep = ""))
        if (u > hmax) 
            y <- (hmax - hdots) * y/u
        y <- hdots + y
    }
    else {
        y <- rep(hdots, length(x))
    }
    if (!is.null(labels)) 
        text(x, y, labels = labels, cex = cex)
    else points(x, y, pch = pch, cex = cex)
    points.coord <- data.frame(x, y)
    if (axes) {
        segments(xmin - b, hdots - h/4, xmax + b, hdots - h/4)
        x <- pretty(x, n = 3, h = 0.5)
        x <- x[(x > xmin - b) & (x < xmax + b)]
        for (i in seq(x)) segments(x[i], hdots - h/4, x[i], hdots - 
            h/2)
        y <- rep(hdots - h, length(x))
        text(x, y, labels = x, cex = cex.axis)
    }
    if (base && !axes) 
        segments(xmin - b, hdots - h/4, xmax + b, hdots - h/4)
    invisible(points.coord)
}
