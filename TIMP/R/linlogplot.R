"linlogplot" <- function (x, y, mu, alpha, xlim, ...) 
{
    maxorigx <- max(x)
    minorigx <- min(x)
    for (i in 1:length(x)) {
        if ((x[i] - mu)/alpha >= 1) 
            x[i] <- alpha + (alpha * log10((x[i] - mu)/alpha))
        else if ((x[i] - mu)/alpha < -1) 
            x[i] <- -alpha - (alpha * log10((mu - x[i])/alpha))
        else x[i] <- x[i] - mu
    }
    for (i in 1:length(xlim)) {
        if ((xlim[i] - mu)/alpha >= 1) 
            xlim[i] <- alpha + (alpha * log10((xlim[i] - mu)/alpha))
        else if ((xlim[i] - mu)/alpha < -1) 
            xlim[i] <- -alpha - (alpha * log10((mu - xlim[i])/alpha))
        else xlim[i] <- xlim[i] - mu
    }

    ticsl <- c(-alpha)
    tics <- c(-alpha)
    cntmin <- -alpha
    while (cntmin > minorigx) {
        cntmin <- cntmin * 10
        ticsl <- append(ticsl, cntmin)
        tics <- append(tics, -alpha - (alpha * log10(-cntmin/alpha)))
    }
    ticsl <- append(sort(ticsl), c(0, alpha))
    tics <- append(sort(tics), c(0, alpha))
    cntmax <- alpha
    while (cntmax < maxorigx) {
        cntmax <- cntmax * 10
        ticsl <- append(ticsl, cntmax)
        tics <- append(tics, alpha + (alpha * log10(cntmax/alpha)))
    }
    plot(x, y, xaxt = "n", xlim = xlim, ...)
    axis(side = 1, at = tics, labels = ticsl)
}
