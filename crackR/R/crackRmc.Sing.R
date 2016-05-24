crackRmc.Sing <-
function(parameters)
{

    Np <- parameters$Np

    ## generate Kc
    kc <- parameters$dta$kc.rsamp(Np)

    ismc.bool <- parameters$ismc.bool

    ## generate eifs
    if( ismc.bool == FALSE )
        {
            eifs   <- parameters$dta$ifs.ractual(Np)
        } else {
            eifs   <- parameters$dta$ifs.rsamp(Np)
            weight <- parameters$dta$ifs.dactual(eifs) / parameters$dta$ifs.dsamp(eifs)
        }

    ## output all data?
    output.all <- parameters$mc.output.all.data

    if( output.all )
        if(ismc.bool)
            {
                all.results <- data.frame(EIFS=eifs,
                                          weight=weight,
                                          Kc=kc,
                                          first.fail=rep(0,Np))
                
            } else {
                all.results <- data.frame(EIFS=eifs,
                                          Kc=kc,
                                          first.fail=rep(0,Np))
            }
    
    
    inspections <- parameters$inspections

    if( ismc.bool && ( length(inspections$flt.interval) > 1 ) ) 
        {
            inspections <- inspections[1,]
            warning("\nthere cannot be inspections with importance sampling.\nrunning to first inspection only.")
        }

    ## flights for interval and record keeping
    flight.int <- inspections$flt.interval
    flight.cum <- cumsum( flight.int )
    n.int <- length( flight.int )

    ## flights already completed at start of each interval
    flight.cum.prev <- c(0, flight.cum[-n.int])

    ms.loc    <- parameters$ms.gumbel[1]
    ms.scale  <- parameters$ms.gumbel[2]

    cg.cc     <- parameters$cg.cc

    ## record the first failure for each trial, (either failure mode)
    first.fail <- rep(max(flight.cum) + 1,Np)

    ## record number of repairs (add one to appropriate slot for each repair)
    repair.count <- rep(0, n.int)

    for(iii in 1:Np)
        {
            
            ## track whether failure occurs; when it does, exit the loop
            ff.bool <- FALSE
            
            ## first fail flight for kc and cc failure modes, starts at one longer than analysis life
            ff.kc <- max(flight.cum) + 1
            ff.cc <- ff.kc

            ## initialize the crack size at end of last interval
            ## here, it's just eifs
            last.crack <- eifs[iii]
            
            for(jjj in 1:n.int)
                {
                    seq.flight <- 1:flight.int[jjj]

                    ## find first flight number
                    first.flight <- approx(x=parameters$dta$cg$crack,
                                           y=parameters$dta$cg$flight,
                                           xout=last.crack, rule=2)$y
                    
                    ## crack length at the end of each flight
                    seq.crack <- approxExtrap(x=parameters$dta$cg$flight,
                                              y=parameters$dta$cg$crack,
                                              xout=seq.flight+first.flight)$y
                    

                    ## find K/Sigma for each flight
                    seq.ksig <- approxExtrap(x=parameters$dta$geo$crack,
                                             y=parameters$dta$geo$ksig,
                                             xout=seq.crack)$y

                    ## ending crack size in sequence (largest)
                    last.crack <- seq.crack[length(seq.crack)]
                    
                    ## test if cg.cc is breached in any flight, if so, find first
                    if( last.crack > cg.cc )
                        {
                            ff.bool <- TRUE
                            ff.cc   <- min( which( seq.crack > cg.cc ) ) + flight.cum.prev[jjj]
                        }
                    
                    ## generate max stress per flight values for each flight
                    seq.ms <- rgumbel(flight.int[jjj], loc=ms.loc, scale=ms.scale)
                    
                    ## K for each flight
                    k <- seq.ksig * seq.ms
                    
                    ## test if Kc is breached in any flight, if so, find first
                    if( max(k) > kc[iii] )
                        {
                            ff.bool <- TRUE
                            ff.kc   <- min( which( k > kc[iii] ) ) + flight.cum.prev[jjj]
                        }
                    
                    if( ff.bool )
                        {
                            first.fail[iii] <- min(ff.cc, ff.kc)
                            break
                        } else if( ismc.bool == FALSE ){
                            
                            ## conduct inspection
                            pod <- parameters$pod.func[[ inspections[jjj,2] ]]( last.crack )
                            if( runif(1) < pod )
                                {
                                    ## count number of repairs for this inspection
                                    repair.count[jjj] <- repair.count[jjj] + 1
                                    ## perform repair
                                    kc[iii]    <- parameters$dta$kc.rsamp(1)
                                    last.crack <- parameters$dta$rfs.ractual(1)
                                }
                        }
                }
        }

    if( output.all ) all.results$first.fail <- first.fail

    
    ## post-process results
    ## use flt.calc.interval to average results within an interval (if =1, will yield all flight results)

    calc.interval <- parameters$flt.calc.interval

    results <- list()
    results$Np.total <- Np
    results$sfpof <- data.frame()
    results$pcd   <- data.frame()
    class(results) <- c("crackRresults", "list")

    ## need to loop again over the intervals
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
                    if( ismc.bool == FALSE )
                        {
                            surv.bool <- ( first.fail >= flights.per.calc.start[iii] )
                            surv <- sum( surv.bool )
                            fail <- sum( surv.bool & ( first.fail <= flights.per.calc.cum[iii] ) )
                            sfpof.this.int[iii] <- ( fail / surv ) / flights.per.calc[iii]
                        } else {
                            surv.bool <- ( first.fail >= flights.per.calc.start[iii] )
                            surv <- sum( weight[ surv.bool ] )
                            fail <- sum( weight[ surv.bool & ( first.fail <= flights.per.calc.cum[iii] ) ] )
                            sfpof.this.int[iii] <- ( fail / surv ) / flights.per.calc[iii]
                        }
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

            results$sfpof$sfpof[results$sfpof$sfpof < parameters$sfpof.min] <- parameters$sfpof.min
            
            ## PCD results
            pcd.temp <- repair.count[jjj] / sum( first.fail > flight.cum[[jjj]] )
            results$pcd <- rbind(results$pcd,
                                 data.frame(flight=flight.cum[[jjj]],
                                            insp.type=inspections$type[jjj],
                                            pcd=pcd.temp))
            
        }

    ## counts for output - in case you want to refine results later
    n.flights <- max(flight.cum) + 1
    fail.count <- rep(0, n.flights)
    for(iii in 1:n.flights)
        fail.count[iii] <- sum( first.fail == iii )
    
    MC.results <- list()
    MC.results$fail.count <- fail.count
    MC.results$repair.count <- repair.count

    if( output.all ) MC.results$all.results <- all.results

    results$MC.results <- MC.results

    return(results)
}
