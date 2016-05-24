calcSfpofFromPofInt <-
function(pof.int)
{
    ## function to format pof.int as SFPOF interval output
    ## used by crackRresults methods if plot parameter sfpof.int = TRUE
    
    flights.int <- pof.int$flights.int
    flights <- cumsum(flights.int)
    n.flights <- length(flights)

    pof.int <- pof.int$pof.int

    ## calculate typical failure prob for a flight in the interval
    ##   note survival to each flight in the interval has already been modeled
    sfpof.int <- 1 - ( ( 1 - pof.int ) ^ (1 / flights.int) )
    
    ## re-format so there's an entry at first and last flight in each interval
    flights.piecewise <- sort( c( c(1,flights[-n.flights]+1) , flights ) )
    sfpof.piecewise <- rep(sfpof.int, each=2)
    sfpof.out <- data.frame(flight=flights.piecewise, sfpof=sfpof.piecewise)
    
    return(sfpof.out)
}
