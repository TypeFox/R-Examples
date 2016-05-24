moment <- function( k=1,           # which moment?
            x,                     # data
            centered=TRUE,         # centered on mean?
            na.rm =T)              # remove missing vals?
            {      

    if (na.rm) { x <- x[!is.na(x)] }

    if (length(k) > 1) {  # vectorize this (fancy)
        return( sapply(k, moment, x=x, centered=centered) )
    }

    if ( centered ) { m = mean(x) } else { m = 0 }

    return ( sum( (x-m)^k ) / length(x) )
}
x <- (1:10)^2; n <- length(x)
moment(1:2,x,centered=F)
moment(1:2,x,centered=T)

c(mean(x), (n-1) / n * var(x))
