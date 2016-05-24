crackRcalcSfpofMc <-
function(fail.count, insp.intervals, calc.interval, sfpof.min=1e-16)
{
    ## function to calculate SFPOF from the raw output of the MC routine
    ## the MC routines should each call this, but at the moment they are redundant
    ## plan to update that in a later version
  
  ## flights for interval and record keeping
  flight.int <- insp.intervals
  flight.cum <- cumsum( flight.int )
  n.int <- length( flight.int )
  
  n.flights <- sum( flight.int ) + 1
  
  ## flights already completed at start of each interval
  flight.cum.prev <- c(0, flight.cum[-n.int])
  
  results <- list()
  class(results) <- c("crackRresults", "list")
  
  results$Np.total <- sum( fail.count )
  
  results$sfpof <- data.frame()
  
  for(jjj in 1:n.int)
  {
      
    ## determine the total number of subintervals for this inspection interval,
    ##   and find the number of flights in each subinterval
    if( flight.int[jjj] < calc.interval )
    {
      n.subint <- 1
    } else {
      n.subint <- round( flight.int[jjj] / calc.interval )
    }
    flights.per.calc.cum <- round(seq(from=flight.int[jjj] / n.subint, to=flight.int[jjj], length=n.subint))
    flights.per.calc <- c(flights.per.calc.cum[1],
                          flights.per.calc.cum[-1] - flights.per.calc.cum[-n.subint])
    
    flights.per.calc.start <- c(1,1+flights.per.calc.cum[-n.subint])
    
    ## update .cum and .start to reflect the actual flight numbers for this interval
    flights.per.calc.cum   <- flights.per.calc.cum   + flight.cum.prev[jjj]
    flights.per.calc.start <- flights.per.calc.start + flight.cum.prev[jjj]
    
    sfpof.this.int <- rep(0, n.subint)
    for(iii in 1:n.subint)
    {
      surv <- sum( fail.count[ flights.per.calc.start[iii] : n.flights ] )
      fail <- sum( fail.count[ flights.per.calc.start[iii] : flights.per.calc.cum[iii] ] )
      sfpof.this.int[iii] <- ( fail / surv ) / flights.per.calc[iii]
    }
    
    if( calc.interval == 1 )
    {
      results.temp <- data.frame(flight=flights.per.calc.cum, sfpof=sfpof.this.int)
    } else {
      results.temp <- data.frame(flight=c(flights.per.calc.start,flights.per.calc.cum),
                                 sfpof=c(sfpof.this.int,sfpof.this.int))
      results.temp <- results.temp[order(results.temp$flight),]
    }
    row.names(results.temp) <- NULL
    
    results$sfpof <- rbind(results$sfpof, results.temp)
    
  }
  
  results$sfpof[results$sfpof < sfpof.min] <- sfpof.min
  
  ## XXX add pof.int to results?
  
  return(results)
}
