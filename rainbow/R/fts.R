`fts` <- function (x, y, start = 1, frequency = 1, xname, yname) 
{
    if (missing(xname)) 
        xname <- deparse(substitute(x))
    if (missing(yname)) 
        yname <- deparse(substitute(y))
    y <- as.matrix(y)
    if (length(x) != nrow(y)) 
        stop("Dimensions do not match")
    ytimes <- time(ts(rep(NA, ncol(y)), start = start, frequency = frequency))
    if(any(colnames(y) == "")){        
        if (max(abs(ytimes - floor(ytimes))) < 1e-09) 
            ylab <- paste(ytimes)
        else if (frequency == 4) 
            ylab <- paste(floor(ytimes), "-Q", cycle(ytimes), sep = "")
        else if (frequency == 12) {
            Mth <- month.abb[cycle(ytimes)]
            ylab <- paste(floor(ytimes), "-", Mth, sep = "")
        }
        else if (frequency == 7) {
            Day <- c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")[cycle(ytimes)]
            ylab <- paste(floor(ytimes), "-", Day, sep = "")
        }
        else ylab <- paste(floor(ytimes), "-", cycle(ytimes), sep = "")
        colnames(y) <- ylab
    }
    else{
        colnames(y) = colnames(y)
    }
    rownames(y) <- x
    return(structure(list(x = x, y = y, time = ytimes, xname = xname, 
        yname = yname), class = c("fts","fds")))

}



