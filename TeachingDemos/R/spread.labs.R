spread.labs <- function(x, mindiff, maxiter=1000, stepsize=1/10,
                          min=-Inf, max=Inf) {
    unsort <- order(order(x))
    x <- sort(x)
    df <- x[-1] - x[ -length(x) ]

    stp <- mindiff * stepsize

    i <- 1
    while( any( df < mindiff ) ) {
        tmp <- c( df < mindiff, FALSE )
        if( tmp[1] && (x[1] - stp) < min ) {  # don't move bottom set
            tmp2 <- as.logical( cumprod(tmp) )
            tmp <- tmp & !tmp2
        }
        x[ tmp ] <- x[ tmp ] - stp
        tmp <- c( FALSE, df < mindiff )
        if( tmp[length(tmp)] && (x[length(x)] + stp) > max ) { # don't move top
            tmp2 <- rev( as.logical( cumprod( rev(tmp) ) ) )
            tmp <- tmp & !tmp2
        }
        x[ tmp ] <- x[ tmp] + stp

        df <- x[-1] - x[-length(x)]
        i <- i + 1
        if( i > maxiter ) {
            warning("Maximum iterations reached")
            break
        }
    }
    x[unsort]
}
