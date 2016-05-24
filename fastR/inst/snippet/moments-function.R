moment <- function(
    k=1,                                    # which moment
    vals=1:6,                               # dice by default
    probs=rep(1/length(vals),length(vals)), # uniform probs
    centered=FALSE) {                       # center on mean of data?

    if (length(k) > 1) {  # vectorize this (fancy)
        return( sapply(k, moment, vals=vals, probs=probs, centered=centered) )
    }

    if ( centered ) {
        m = sum(vals * probs)
    } else { 
        m = 0
    }
    sum((vals-m)^k * probs)
}

moment(k=1,0:10,dbinom(0:10,10,0.4))
moment(k=2,0:10,dbinom(0:10,10,0.4), centered=FALSE)
moment(k=2,0:10,dbinom(0:10,10,0.4), centered=TRUE)
10 * 0.4 * 0.6   # should match previous and next value
moment(k=2,0:10,dbinom(0:10,10,0.4),centered=FALSE) - 
    moment(k=1, 0:10,dbinom(0:10,10,0.4), centered=FALSE)^2
round(moment(k=1:4, 0:10,dbinom(0:10,10,0.4), centered=FALSE),5)
round(moment(k=1:4, 0:10,dbinom(0:10,10,0.4), centered=TRUE),5)
