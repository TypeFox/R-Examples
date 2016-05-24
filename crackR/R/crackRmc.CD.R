crackRmc.CD <-
function(parameters)
{
  
  Np <- parameters$Np
  
  ## generate Kc (use primary cold since all are same)
  kc <- parameters$dta.pc$kc.rsamp(Np)

  ## whether this is an importance sampling run
  ismc.bool <- parameters$ismc.bool

  ## for partitioning PCD (new feature)
  pod.threshold        <- parameters$pod.threshold
  n.rep.types          <- length(pod.threshold) + 1
  two.hole.repair.both <- parameters$two.hole.repair.both
  two.hole.repair.both <- c(FALSE, two.hole.repair.both) ## never repair both for smallest type repair
  
  ## generate eifs for primary and secondary cracks, independent
  if( ismc.bool == FALSE )
  {
    eifs.p   <- parameters$dta.pc$ifs.ractual(Np)
    eifs.s   <- parameters$dta.sc$ifs.ractual(Np)  
  } else {
    eifs.p   <- parameters$dta.pc$ifs.rsamp(Np)
    eifs.s   <- parameters$dta.sc$ifs.rsamp(Np)
    weight <- ( parameters$dta.pc$ifs.dactual(eifs.p) * parameters$dta.sc$ifs.dactual(eifs.s) ) / 
      ( parameters$dta.pc$ifs.dsamp(eifs.p) * parameters$dta.sc$ifs.dsamp(eifs.s) )
  }

    ## output all data?
    output.all <- parameters$mc.output.all.data

    if( output.all )
        if(ismc.bool)
            {
                all.results <- data.frame(EIFS.p=eifs.p,
                                          EIFS.s=eifs.s,
                                          weight=weight,
                                          Kc=kc,
                                          first.fail=rep(0,Np))
                
            } else {
                all.results <- data.frame(EIFS.p=eifs.p,
                                          EIFS.s=eifs.s,
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
  
  ms.loc   <- parameters$ms.gumbel[1]
  ms.scale <- parameters$ms.gumbel[2]
  
  cg.cc.pc <- parameters$cg.cc.pc
  cg.cc.ph <- parameters$cg.cc.ph
  cg.cc.sc <- parameters$cg.cc.sc
  cg.cc.sh <- parameters$cg.cc.sh
  
  ## record the first failure for each trial, (either failure mode, either crack)
  first.fail <- rep(max(flight.cum) + 1,Np)
  
  ## record number of repairs (add one to appropriate slot for each repair)
  ## for an importance sampling run these will stay zeros
  repair.count.p <- rep(0, n.int)
  repair.count.s <- rep(0, n.int)
  
  repair.count.by.type <- data.frame( matrix(0, nrow=n.int, ncol=n.rep.types) )
  names(repair.count.by.type) <- paste("pcd", 1:n.rep.types, sep="")
  
  for(iii in 1:Np)
  {
      
      ## track whether failure occurs; when it does, exit the loop
      ff.bool <- FALSE
      
      ## first fail flight for kc and cc failure modes, starts at one longer than analysis life
      ff.kc <- max(flight.cum) + 1
      ff.cc <- ff.kc

      ## initialize the crack size at end of last interval
      ## before we begin, it's just eifs
      last.crack.p <- eifs.p[iii]
      last.crack.s <- eifs.s[iii]
      
      ## boolean which identifies whether each crack is hot yet (always cold to start)
      ## we could use two for primary and secondary, but currently just using one
      ##     hot.bool.p <- FALSE
      ##     hot.bool.s <- FALSE
      hot.bool <- FALSE
      
      for(jjj in 1:n.int)
          {
              seq.flight <- 1:flight.int[jjj]
              
              ## need to grow each crack until the other goes critical, then switch curves
              
              ## possible way to do this:
              ##   if either pri. or sec. is already hot, just use hot tables
              ##   if both cold, use cold tables to find if and where a crack will get hot
              ##     if neither will get hot in the int, use cold tables
              ##     if one will get hot, run with cold tables until that flight, then hot after
              
              if( hot.bool )
                  {

                      ## find first flight number (hot tables)
                      first.flight.p <- approx(x=parameters$dta.ph$cg$crack,
                                               y=parameters$dta.ph$cg$flight,
                                               xout=last.crack.p, rule=2)$y
                      first.flight.s <- approx(x=parameters$dta.sh$cg$crack,
                                               y=parameters$dta.sh$cg$flight,
                                               xout=last.crack.s, rule=2)$y

                      ## either crack is hot (this could be faster by only growing the hot crack...)
                      seq.crack.p <- approxExtrap(x=parameters$dta.ph$cg$flight,
                                                  y=parameters$dta.ph$cg$crack,
                                                  xout=seq.flight+first.flight.p)$y
                      seq.crack.s <- approxExtrap(x=parameters$dta.sh$cg$flight,
                                                  y=parameters$dta.sh$cg$crack,
                                                  xout=seq.flight+first.flight.s)$y
                      
                      seq.ksig.p <- approx(x=parameters$dta.ph$geo$crack,
                                           y=parameters$dta.ph$geo$ksig,
                                           xout=seq.crack.p, rule=2)$y
                      seq.ksig.s <- approx(x=parameters$dta.sh$geo$crack,
                                           y=parameters$dta.sh$geo$ksig,
                                           xout=seq.crack.s, rule=2)$y
                      
                  } else {
                      ## neither is hot, so grow using cold, then check if one becomes hot in int

                      ## find first flight number (cold tables)
                      first.flight.p <- approx(x=parameters$dta.pc$cg$crack,
                                               y=parameters$dta.pc$cg$flight,
                                               xout=last.crack.p, rule=2)$y
                      first.flight.s <- approx(x=parameters$dta.sc$cg$crack,
                                               y=parameters$dta.sc$cg$flight,
                                               xout=last.crack.s, rule=2)$y
                      
                      seq.crack.p <- approxExtrap(x=parameters$dta.pc$cg$flight,
                                                  y=parameters$dta.pc$cg$crack,
                                                  xout=seq.flight+first.flight.p)$y
                      seq.crack.s <- approxExtrap(x=parameters$dta.sc$cg$flight,
                                                  y=parameters$dta.sc$cg$crack,
                                                  xout=seq.flight+first.flight.s)$y
                      
                      seq.ksig.p <- approx(x=parameters$dta.pc$geo$crack,
                                           y=parameters$dta.pc$geo$ksig,
                                           xout=seq.crack.p, rule=2)$y
                      seq.ksig.s <- approx(x=parameters$dta.sc$geo$crack,
                                           y=parameters$dta.sc$geo$ksig,
                                           xout=seq.crack.s, rule=2)$y
                      
                      if( max(seq.crack.p) >= cg.cc.pc || max(seq.crack.s) >= cg.cc.sc )
                          {
                              ## at least one of the cracks reaches critical in the interval
                              hot.bool <- TRUE
                              
                              ## identify which flights occur before one goes hot
                              cold.flights <- ( seq.crack.p <= cg.cc.pc ) & ( seq.crack.s <= cg.cc.sc )

                              ## if it goes hot in the very first flight, the below will error
                              ##   force the first flight to be cold
                              if( max(cold.flights) == 0 ) ## all are false, i.e., all flights are hot
                                  cold.flights[1] <- TRUE
                              
                              ## keep the crack sizes and K/sigma values for these cold flights
                              seq.crack.p.1 <- seq.crack.p[ cold.flights ]
                              seq.crack.s.1 <- seq.crack.s[ cold.flights ]
                              
                              seq.ksig.p.1 <- seq.ksig.p[ cold.flights ]
                              seq.ksig.s.1 <- seq.ksig.s[ cold.flights ]
                              
                              ## need to update first.flight.p/s to correspond to hot tables
                              first.flight.p <- approx(x=parameters$dta.ph$cg$crack,
                                                       y=parameters$dta.ph$cg$flight,
                                                       xout=seq.crack.p.1[length(seq.crack.p.1)],
                                                       rule=2)$y
                              first.flight.s <- approx(x=parameters$dta.sh$cg$crack,
                                                       y=parameters$dta.sh$cg$flight,
                                                       xout=seq.crack.s.1[length(seq.crack.s.1)],
                                                       rule=2)$y
                              
                              ## this shorter flight sequence will be used with hot tables
                              seq.flight  <- 1:sum(!cold.flights)
                              
                              ## use hot tables to grow both cracks using hot data
                              ##   could possibly speed this by only growing the hot one
                              seq.crack.p.2 <- approxExtrap(x=parameters$dta.ph$cg$flight,
                                                            y=parameters$dta.ph$cg$crack,
                                                            xout=seq.flight+first.flight.p)$y
                              seq.crack.s.2 <- approxExtrap(x=parameters$dta.sh$cg$flight,
                                                            y=parameters$dta.sh$cg$crack,
                                                            xout=seq.flight+first.flight.s)$y
                              
                              seq.ksig.p.2 <- approx(x=parameters$dta.ph$geo$crack,
                                                     y=parameters$dta.ph$geo$ksig,
                                                     xout=seq.crack.p.2, rule=2)$y
                              seq.ksig.s.2 <- approx(x=parameters$dta.sh$geo$crack,
                                                     y=parameters$dta.sh$geo$ksig,
                                                     xout=seq.crack.s.2, rule=2)$y
                              
                              ## combine cold and hot flights
                              seq.crack.p <- c( seq.crack.p.1, seq.crack.p.2 )
                              seq.crack.s <- c( seq.crack.s.1, seq.crack.s.2 )
                              seq.ksig.p  <- c( seq.ksig.p.1,  seq.ksig.p.2 )
                              seq.ksig.s  <- c( seq.ksig.s.1,  seq.ksig.s.2 )
                          }
                  }
              
              ## ending crack size in sequence (largest)
              last.crack.p <- seq.crack.p[ length(seq.crack.p) ]
              last.crack.s <- seq.crack.s[ length(seq.crack.s) ]
              
              ## test if cg.cc is breached BY BOTH CRACKS
              ##   if both breached, find flight when both were
              if( ( last.crack.p > cg.cc.ph ) && ( last.crack.s > cg.cc.sh ) )
                  {
                      ff.bool <- TRUE
                      ff.cc.p <- min( which( seq.crack.p > cg.cc.ph ) )
                      ff.cc.s <- min( which( seq.crack.s > cg.cc.sh ) )
                      ff.cc   <- max( ff.cc.p, ff.cc.s ) + flight.cum.prev[jjj]
                  }
              
              ## generate max stress per flight values for each flight
              seq.ms <- rgumbel(flight.int[jjj], loc=ms.loc, scale=ms.scale)
              
              ## K for each flight
              ## for K/sigma, both cracks need to fail, so only need LOWER value
              seq.ksig.low <- apply( cbind(seq.ksig.p,seq.ksig.s), 1, min )
              k <- seq.ksig.low * seq.ms
              
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
                      pod.p <- parameters$pod.func[[ inspections[jjj,2] ]]( last.crack.p )
                      pod.s <- parameters$pod.func[[ inspections[jjj,2] ]]( last.crack.s )
                      rep.p <- ( runif(1) < pod.p )
                      rep.s <- ( runif(1) < pod.s )

                      ## repair type that would occur if found
                      rep.type.p <- findInterval(last.crack.p, pod.threshold) + 1
                      rep.type.s <- findInterval(last.crack.s, pod.threshold) + 1

                      ## if a repair is performed that fixed both cracks, the more severe occurs
                      rep.type.larger <- max(rep.type.p, rep.type.s)
                      
                      if( rep.p || rep.s )
                          {
                              ## some sort of repair occurs, but are there one or two repairs to do?
                              
                              if( parameters$two.hole.bool )
                                  {
                                      ## there are three variables we need here (using .p)
                                      ##   rep.p : boolean, was the primary crack found?
                                      ##   rep.type.p : int, what type of repair for primary?
                                      ##   rep.type.p.both : boolean, is primary repair type a both fixer? (new)
                                      ##   repair.count.p[jjj] : count number of repairs to primary crack

                                      rep.type.p.both <- two.hole.repair.both[ rep.type.p ]
                                      rep.type.s.both <- two.hole.repair.both[ rep.type.s ]

                                      ## treat these 3 cases separately: p.only, s.only, both

                                      ## p.only
                                      if( rep.p && !rep.s )
                                      {
                                          ## count repairs by crack
                                          repair.count.p[jjj] <- repair.count.p[jjj] + 1
                                          
                                          last.crack.p <- parameters$dta.pc$rfs.ractual(1)

                                          ## if this is a both fixer, only count the one repair,
                                          ##   but fix both cracks -- use .larger in case other is bigger!
                                          if( rep.type.p.both )
                                          {
                                              ## count repairs by repair type (use larger!)
                                              repair.count.by.type[jjj, rep.type.larger] <-
                                                  repair.count.by.type[jjj, rep.type.larger] + 1
                                              ## repair secondary also
                                              last.crack.s <- parameters$dta.sc$rfs.ractual(1)  
                                          } else {
                                              ## count repairs by repair type
                                              repair.count.by.type[jjj, rep.type.p] <-
                                                  repair.count.by.type[jjj, rep.type.p] + 1
                                          }
                                      }
                                      
                                      ## s.only
                                      if( !rep.p && rep.s )
                                      {
                                          ## count repairs by crack
                                          repair.count.s[jjj] <- repair.count.s[jjj] + 1
                                          last.crack.s <- parameters$dta.sc$rfs.ractual(1)

                                          ## if this is a both fixer, only count the one repair,
                                          ##   but fix both cracks -- use .larger in case other is bigger!
                                          if( rep.type.s.both )
                                          {
                                              ## count repairs by repair type (use larger!)
                                              repair.count.by.type[jjj, rep.type.larger] <-
                                                  repair.count.by.type[jjj, rep.type.larger] + 1
                                              ## repair primary also
                                              last.crack.p <- parameters$dta.pc$rfs.ractual(1)
                                          } else {
                                              ## count repairs by repair type
                                              repair.count.by.type[jjj, rep.type.s] <-
                                                  repair.count.by.type[jjj, rep.type.s] + 1
                                          }
                                      }
                                      
                                      ## p and s were found
                                      if( rep.p && rep.s )
                                      {
                                          last.crack.p <- parameters$dta.pc$rfs.ractual(1)
                                          last.crack.s <- parameters$dta.sc$rfs.ractual(1)

                                          kc <- parameters$dta.pc$kc.rsamp(1)

                                          ## if this is a both fixer FOR EITHER, only count one repair
                                          if( rep.type.p.both || rep.type.s.both )
                                          {
                                              ## count repairs by repair type (use larger!)
                                              repair.count.by.type[jjj, rep.type.larger] <-
                                                  repair.count.by.type[jjj, rep.type.larger] + 1
                                              
                                              ## count repairs by crack (USE HALF A REPAIR)
                                              ## the purpose of these repair counts is to report the probability
                                              ## of FINDING AND REPAIRING a crack. since both were found here,
                                              ## it seems to make sense to split it up.
                                              repair.count.p[jjj] <- repair.count.p[jjj] + 0.5
                                              repair.count.s[jjj] <- repair.count.s[jjj] + 0.5
                                          } else {
                                              ## count repairs by repair type
                                              repair.count.by.type[jjj, rep.type.p] <-
                                                  repair.count.by.type[jjj, rep.type.p] + 1
                                              repair.count.by.type[jjj, rep.type.s] <-
                                                  repair.count.by.type[jjj, rep.type.s] + 1
                                              
                                              ## count repairs by crack
                                              repair.count.p[jjj] <- repair.count.p[jjj] + 1
                                              repair.count.s[jjj] <- repair.count.s[jjj] + 1
                                          }
                                      }

                                  } else {

                                      ## for one hole model, the larger repair is assumed to occur
                                      repair.count.by.type[jjj, rep.type.larger] <-
                                          repair.count.by.type[jjj, rep.type.larger] + 1

                                      repair.count.p[jjj] <- repair.count.p[jjj] + 1
                                      repair.count.s[jjj] <- repair.count.s[jjj] + 1
                                      last.crack.p <- parameters$dta.pc$rfs.ractual(1)
                                      last.crack.s <- parameters$dta.sc$rfs.ractual(1)

                                      kc <- parameters$dta.pc$kc.rsamp(1)
                                  }

                              ## need to reset hot.bool since a repair was performed
                              if( ( last.crack.p >= cg.cc.pc ) || ( last.crack.s >= cg.cc.sc ) )
                                  {
                                      hot.bool <- TRUE
                                  } else {
                                      hot.bool <- FALSE
                                  }
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
  results$pcd.old   <- data.frame()
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
    if( parameters$two.hole.bool )
    {
      ## allowed to double count repairs since PCD is really the exp. # of rep's here
      ## that is, in the two hole model, PCD can be > 1
      ##   this allows us to use the same routines for calculating expected costs
      pcd.temp <- (repair.count.p[jjj] + repair.count.s[jjj]) / sum( first.fail > flight.cum[[jjj]] )
    } else {
      ## repair.count.p equals repair.count.s for one hole model; only need one
      pcd.temp <- repair.count.p[jjj] / sum( first.fail > flight.cum[[jjj]] )    
    }
    results$pcd.old <- rbind(results$pcd.old,
                         data.frame(flight=flight.cum[[jjj]],
                                    insp.type=inspections$type[jjj],
                                    pcd=pcd.temp))

    pcd.temp <- repair.count.by.type[jjj,] / sum( first.fail > flight.cum[[jjj]] )
    names(pcd.temp) <- paste("pcd", 1:n.rep.types, sep="")
    results$pcd <- rbind(results$pcd,
                         data.frame(flight=flight.cum[[jjj]],
                                    insp.type=inspections$type[jjj],
                                    pcd.temp,
                                    pcd=sum(pcd.temp)))    
  }
  
  ## create a crackRmc object so we can add runs to an existing run?
  ## for now just return the results

  ## counts for output - in case you want to refine results later
  n.flights <- max(flight.cum) + 1
  fail.count <- rep(0, n.flights)
  for(iii in 1:n.flights)
      fail.count[iii] <- sum( first.fail == iii )
      
  MC.results <- list()
  MC.results$fail.count <- fail.count
  MC.results$repair.count.p <- repair.count.p
  MC.results$repair.count.s <- repair.count.s

  if( output.all ) MC.results$all.results <- all.results
  
  return(results)
}
