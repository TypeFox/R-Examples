disana <- function (x,panel='all') 
{
    if (class(x) == "dist") {
        y <- as.matrix(x)
        triang <- x
    } else {
        y <- as.matrix(x)
        triang <- y[row(y) > col(y)]
    }
    is.na(diag(y)) <- TRUE
    tmin <- apply(y, 1, function(z) {
        min(z, na.rm = TRUE)
    })
    tavg <- apply(y, 1, function(z) {
        mean(z, na.rm = TRUE)
    })
    tmax <- apply(y, 1, function(z) {
        max(z, na.rm = TRUE)
    })
    plots <- NULL
    if (panel=='all' || panel==1) {
        plot(sort(triang), xlab = "Sorted Value", ylab = "Dissimilarity")
        if (panel=='all') readline("Press return for next page....")
    }
    if (panel=='all' || panel==2) {
        plot(sort(tmin), ylim = c(0, max(tmax)), xlab = "Sorted Plot", 
            ylab = "Dissimilarity")
        points(sort(tavg), col = 2)
        points(sort(tmax), col = 3)
        if (panel=='all') readline("Press return for next page....")
    }
    if (panel=='all' || panel==3) {
        plot(tmin, tavg, xlab = "Minimum Dissimilarity", ylab = "Average Dissimilarity")
        lines(c(0.5, 0.5), c(min(tavg), max(tavg)), col = 2)
        yorn <- readline("Do you want to identify individual plots [Y or N] : ")
        if (yorn == "Y" || yorn == "y") 
            plots <- identify(tmin, tavg, attr(x, "Labels"))
    }
    res <- list(min = tmin, mean = tavg, max = tmax, plots = plots)
    attr(res,'call') <- match.call()
    invisible(res)
}

