`getCentroid` <-
function(data) {
    return(lapply(data, function(x) {
        if (all(is.na(x[nrow(x), ]))) {
            x <- x[-nrow(x), ]
        }
        
        if (nrow(x) > 10) {
            vals <- seq(1, nrow(x), by = floor(nrow(x)/10))
            x <- x[c(vals, nrow(x)), ]
        }
        
        asum <- cxsum <- cysum <- 0
        
        for (i in seq(2, nrow(x))) {
            asum <- asum + x$x[i - 1]*x$y[i] - x$x[i]*x$y[i - 1]
            cxsum <- cxsum + (x$x[i - 1] + x$x[i])*(x$x[i - 1]*x$y[i] - x$x[i]*x$y[i - 1])
            cysum <- cysum + (x$y[i - 1] + x$y[i])*(x$x[i - 1]*x$y[i] - x$x[i]*x$y[i - 1])
        }
        
        return(c((1/(3*asum))*cxsum, (1/(3*asum))*cysum))
    }))
}


