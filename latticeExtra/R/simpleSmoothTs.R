

simpleSmoothTs <-
    function(x, ...)
{
    UseMethod("simpleSmoothTs")
}

simpleSmoothTs.default <-
    function(x, ...,
             width = NROW(x) %/% 10 + 1, n = NROW(x),
             c = 1, sides = 2, circular = FALSE,
             kern = kernel("daniell", rep(floor((width/sides)/sqrt(c)), c)))
{

    if (sides == 2) {
        ii <- -kern$m:kern$m
        filter <- kern[ii]
    } else if (sides == 1) {
        ii <- 0:kern$m
        filter <- kern[ii] / sum(kern[ii]) ## normalise
    } else stop("unrecognised value of 'sides'")
    x <- as.ts(x)
    xf <- x
    xf[] <- filter(x, filter, sides = sides, circular = circular)
    if (n < NROW(x)) {
        ## reduce the number of points by aggregating chunks of 'reduce' time steps
        reduce <- round(NROW(x) / n)
        if (reduce > 1) {
            ndeltat <- deltat(xf) * reduce
            ## work-around for bug in aggregate.ts
            if ((ndeltat > 1) && (getRversion() < "2.11.1"))
                ndeltat <- ndeltat * (1 + getOption("ts.eps")/1000)
            xf <- aggregate(xf, ndeltat = ndeltat, FUN = mean)
            ## and adjust it so that each point is centered compared to the original series
            tsp(xf)[1:2] <- tsp(xf)[1:2] + (deltat(xf) %/% 2)
        }
    }
    xf
}

simpleSmoothTs.zoo <-
    function(x, ..., n = NROW(x))
{
    xts <- as.ts(x)
    xtsf <- simpleSmoothTs(xts, ..., n = n)
    ii <- TRUE
    if (n < NROW(x)) {
        ## find indices of aggregated time series in the original data
        ii <- findInterval(time(xtsf), time(xts))
        ## extract elements of time index corresponding to aggregated series
        ans <- zoo::zoo(as.matrix(xtsf), time(x)[ii])
    } else {
        ans <- x
        zoo::coredata(ans) <- zoo::coredata(xtsf)
    }
    ans
}
