`sfts` <-
function (data, period = frequency(data), start = tsp(data)[1], 
    frequency = 1, xname, yname) 
{
    if (class(data) != "ts") {
        warning("object is not time series.")
    }
    else {
        y = matrix(data, period, )
        x = 1:period
        if (missing(xname)) 
            xname <- deparse(x)
        if (missing(yname)) 
            yname <- "sliced data"
        if (length(x) != nrow(y)) 
            stop("Dimensions do not match")
        ytimes <- time(ts(rep(NA, ncol(y)), start = start, frequency = frequency))
        if (max(abs(ytimes - floor(ytimes))) < 1e-09) 
            ylab <- paste(ytimes)
        else if (frequency == 4) 
            ylab <- paste(floor(ytimes), "-Q", cycle(ytimes), 
                sep = "")
        else if (frequency == 12) {
            Mth <- month.abb[cycle(ytimes)]
            ylab <- paste(floor(ytimes), "-", Mth, sep = "")
        }
        else if (frequency == 7) {
            Day <- c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", 
                "Sat")[cycle(ytimes)]
            ylab <- paste(floor(ytimes), "-", Day, sep = "")
        }
        else ylab <- paste(floor(ytimes), "-", cycle(ytimes), 
            sep = "")
        colnames(y) <- ylab
        rownames(y) <- x
        return(structure(list(x = x, y = y, time = ytimes, xname = xname, 
            yname = yname), class = c("sfts", "fts", "fds")))
    }
}
