moment.cont <- function( k=1,       # which moment?
        dist = dnorm,  
        args=list(),                # arguments to dist()
        range=c(-Inf,Inf),          
        centered=FALSE) {           # centered on mean?

    if (length(k) > 1) {  # vectorize this (fancy)
        return( sapply(k, moment.cont, 
                        dist=dist, args=args,
                        range=range, centered=centered) )
    }

    if ( centered ) {
        m = moment.cont(dist=dist, range=range, k=1, centered=FALSE)
    } else { 
        m = 0
    }
    int.out <- integrate(
                    function(x) { (x-m)^k * dist(x) }, 
                    range[1], range[2])
    return (int.out$value)
}

moment.cont(dunif,k=1,centered=FALSE)
moment.cont(dunif,k=2,centered=FALSE)
moment.cont(dunif,k=2,centered=TRUE)
moment.cont(dunif,k=1:4,centered=FALSE)
round(moment.cont(dunif,k=1:4,centered=TRUE),5)
round(moment.cont(dnorm,k=1:4,centered=TRUE),5)
round( moment.cont(function(x) {dnorm(x,10,3)},k=1:4,centered=TRUE), 5)
