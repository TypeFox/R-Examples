`ARToMA` <-
function(phi, lag.max)
{
    p <- length(phi)
    x <- numeric(lag.max + 1)
    x <- 1
    for(i in 1:p) {
        x[i + 1] <- crossprod(phi[1:i], (rev(x))[1:i])
    }
    if(lag.max > p) {
        for(i in (p + 1):lag.max) {
            x[i + 1] <- crossprod(phi, (rev(x))[1:p])
        }
    }
    return(x)
}

