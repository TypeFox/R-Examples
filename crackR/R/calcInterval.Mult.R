calcInterval.Mult <-
function(obj, interval.flights)
{
    ## takes a crackR Mult object and advances a specified interval of flights
    ## this uses multiple DTA types and is slower than the Sing routine
    
    insp.num      <- dim(obj$results$pcd)[1] + 1      ## current insp number
    calc.interval <- obj$parameters$flt.calc.interval ## how often to output results

    ## determine the total number of subintervals for this inspection interval,
    ##   and find the number of flights in each subinterval
    if( interval.flights < calc.interval )
        {
            n.subint <- 1
        } else {
            n.subint <- round( interval.flights / calc.interval )
        }
    flights.per.calc.cum <- round(seq(from=interval.flights / n.subint, to=interval.flights, length=n.subint))
    flights.per.calc <- c(flights.per.calc.cum[1],
                          flights.per.calc.cum[-1] - flights.per.calc.cum[-n.subint])
    ## each component of flights.per.calc is the number of flights for that subinterval calculation

    ## number of flights already contained in obj$results
    flights.completed <- ifelse(is.null(obj$results$sfpof$flight), 0, max(obj$results$sfpof$flight))

    ## identify the actual flight numbers at which results will be appended
    calc.flights <- flights.completed +
        as.numeric(rbind(flights.per.calc.cum + 1 - flights.per.calc, flights.per.calc.cum))

    dta       <- obj$parameters$dta       ## deterministic damage tolerance analysis data
    dta.types <- obj$parameters$dta.types ## how many distinct sets are in the dta block?

    ms.lc     <- obj$parameters$ms.gumbel[1] ## gumbel max stress dist is hard coded at present
    ms.sc     <- obj$parameters$ms.gumbel[2] ## " "
    sfpof.min <- obj$parameters$sfpof.min    ## smaller estimates of sfpof are truncated at this value
    cg.cc     <- obj$parameters$cg.cc        ## critical crack length at which pr(fail)=100%
    
    calc.num    <- length(flights.per.calc) ## number of subintervals for this set of calculations
    sfpof.first <- rep(0, calc.num) ## storage of sfpof at start of each subinterval
    sfpof.last  <- rep(0, calc.num) ## storage of sfpof at end " "
    pof.int     <- rep(0, calc.num) ## storage of probability of failure for the subinterval

    ## prepare for bootstrapping of sfpof results
    sfpof.boot.boolean <- obj$parameters$bootstrap.sfpof
    if( sfpof.boot.boolean )
        {
            ## calc.num is number of subintervals, will do bootstrap at the
            ##   start and end of each subinterval (makes better plots this way)
            sfpof.boot <- matrix(0, nrow=2*calc.num, ncol=length(obj$parameters$bootstrap.quantiles))
            colnames(sfpof.boot) = paste("quantile_", obj$parameters$bootstrap.quantiles, sep="")
            sfpof.boot <- as.data.frame(sfpof.boot)
            
            boot.n <- obj$parameters$bootstrap.samples
            boot.q <- obj$parameters$bootstrap.quantiles
        }

    ## model state at the start of the routine
    a   <- obj$state$a
    kc  <- obj$state$kc
    w   <- obj$state$w
    typ <- obj$state$typ
    Np  <- obj$parameters$Np

    ## split calculation here for flight-by-flight or interval calculation

    if( tolower(obj$parameters$survival.updating) == "fbf" )
        {
            ## each particle may be of a different type. identify the types and loop
            f.current    <- rep(0, Np)
            ksig.current <- f.current ## placeholder needed since masking below
            dt.index <- list()
            if(dta.types==1) dt.index[[1]] <- rep(TRUE, Np) else
            for(kkk in 1:dta.types) dt.index[[kkk]] <- ( typ == kkk )

            dta.types.current <- unique(typ) ## these are the types currently present in the state
            for(kkk in dta.types.current)
                f.current[dt.index[[kkk]]] <- approxExtrap(x=dta[[kkk]]$cg$crack,
                                                           y=dta[[kkk]]$cg$flight,
                                                           xout=a[dt.index[[kkk]]])$y

            ## give warning in the loop if any weights are NaN, but only display once; so need flag
            w.NaN <- FALSE

            ## for each calc interval, loop over all flights in the interval and
            ##   save sfpof.first, sfpof.last, a.last, and updated weights
            for(iii in 1:calc.num)
                {
                    ## run a check to see if any weights are NaN since that usually
                    ##   means all particles have reached the critical crack length
                    if(w.NaN == FALSE)
                        if(any(is.na(w)))
                            {
                                w.NaN <- TRUE
                                warning("at least one weight is NaN. this usually means all particles have reached the critical crack length.")
                            }

                    sfpof.temp <- rep(0, flights.per.calc[iii])
                    
                    for(jjj in 1:flights.per.calc[iii])
                        {

                            f.current    <- f.current + 1
                            for(kkk in dta.types.current)
                                {
                                    a[dt.index[[kkk]]] <- approxExtrap(x=dta[[kkk]]$cg$flight,
                                                                       y=dta[[kkk]]$cg$crack,
                                                                       xout=f.current[dt.index[[kkk]]])$y
                                    ksig.current[dt.index[[kkk]]] <- approx(x=dta[[kkk]]$geo$crack,
                                                                            y=dta[[kkk]]$geo$ksig,
                                                                            xout=a[dt.index[[kkk]]],
                                                                            rule=2)$y
                                }

                            ## testing this line; safer than setting them equal to cg.cc (floating point...)
                            a[ a >= cg.cc ] <- Inf

                            ## K>Kc failure mode (DTA type irrelevant)
                            pos.current  <- pgumbel(kc, loc=(ms.lc*ksig.current), scale=(ms.sc*ksig.current)) 
                            pos.current[ a >= cg.cc ] <- 0 ## a>ac failure mode

                            sfpof.temp[jjj] <- 1 - sum(w * pos.current)

                            ## bootstrap for SFPOF at first and last flight in the subinterval
                            if( sfpof.boot.boolean )
                                if( jjj == 1 || jjj == flights.per.calc[iii] )
                                    {
                                        sfpof.boot.n <- rep(0, boot.n)
                                        Np.vec <- 1:Np
                                        for(bbb in 1:boot.n)
                                            {
                                                index.boot <- sample(x=Np.vec, size=Np, replace=TRUE)
                                                pos.boot   <- pos.current[index.boot]
                                                w.boot     <- w[index.boot]
                                                w.boot     <- w.boot / sum(w.boot)

                                                sfpof.boot.n[bbb] <- 1-sum(w.boot*pos.boot)
                                            }
                                        if( jjj == 1)
                                            {
                                                sfpof.boot[2*iii-1,] <- quantile(x=sfpof.boot.n, probs=boot.q)
                                            } else {
                                                sfpof.boot[2*iii,] <- quantile(x=sfpof.boot.n, probs=boot.q)
                                            }
                                    }

                            w <- w * pos.current
                            w <- w / sum(w)   
                        }
                    sfpof.first[iii] <- sfpof.temp[1]
                    sfpof.last[iii]  <- sfpof.temp[flights.per.calc[iii]]
                    pof.int[iii]     <- 1-prod(1-sfpof.temp)
                }

            obj$state$a <- a
            obj$state$w <- w

        } else if( tolower(obj$parameters$survival.updating) == "int" )
            {

                a.frst <- a

                ## identify the particle types
                dt.index <- list()
                if(dta.types==1) dt.index[[1]] <- rep(TRUE, Np) else
                for(kkk in 1:dta.types) dt.index[[kkk]] <- ( typ == kkk )

                ## initialize the various vectors that are filled with masking in loop
                flt.frst <- rep(0, Np)
                flt.mddl <- flt.frst
                flt.last <- flt.frst
                a.mddl <- flt.frst
                a.last <- flt.frst
                ksig.frst <- flt.frst
                ksig.mddl <- flt.frst
                ksig.last <- flt.frst
                a.end.of.last <- flt.frst
                
                dta.types.current <- unique(typ) ## these are the types that are currently present in the state
                for(kkk in dta.types.current)
                    flt.frst[dt.index[[kkk]]] <- approxExtrap(x=dta[[kkk]]$cg$crack,
                                                              y=dta[[kkk]]$cg$flight,
                                                              xout=a.frst[dt.index[[kkk]]])$y

                ## find the flight number for cg.cc; used to update flight hour for cc failed cracks
                ## cg.cc.flight <- rep(0, Np)
                ## for(kkk in dta.types.current)
                ##     cg.cc.flight[dt.index[[kkk]]] <- approxExtrap(x=dta[[kkk]]$cg$crack,
                ##                                                   y=dta[[kkk]]$cg$flight,
                ##                                                   xout=cg.cc)$y

                ## giving a warning in the loop if any weights are NaN
                ## this is a flag since i only want the warning once
                w.NaN <- FALSE
                
                ## for each calc interval, loop over all flights in the interval and save
                ##   sfpof.frst, sfpof.last, a.last, and updated importance weights
                for(iii in 1:calc.num)
                    {

                        ## run a check to see if any weights are NaN since that usually
                        ##   means all particles have reached the critical crack length
                        if(w.NaN == FALSE)
                            if(any(is.na(w)))
                                {
                                    w.NaN <- TRUE
                                    warning("at least one weight is NaN. this usually means all particles have reached the critical crack length.")
                                }
                        
                        ## frst = first flight in the interval
                        ## mddl = middle flight in int
                        ## last = last flight in int

                        num.flights <- flights.per.calc[iii]

                        ##-------------------##
                        ## K>Kc failure mode ##
                        ##-------------------##
                        
                        ## flt.frst already obtained
                        flt.mddl <- flt.frst + num.flights*0.5 + 0.5 ## middle of middle flight
                        flt.last <- flt.frst + num.flights + 1       ## end of last flight

                        ## loop over dta.types
                        for(kkk in dta.types.current)
                            {
                                temp.index <- dt.index[[kkk]]
                                
                                temp.cg.flight <- dta[[kkk]]$cg$flight
                                temp.cg.crack  <- dta[[kkk]]$cg$crack
                                temp.geo.crack <- dta[[kkk]]$geo$crack
                                temp.geo.ksig  <- dta[[kkk]]$geo$ksig
                                
                                a.mddl[temp.index] <- approxExtrap(x=temp.cg.flight,
                                                                   y=temp.cg.crack,
                                                                   xout=flt.mddl[temp.index])$y
                                a.last[temp.index] <- approxExtrap(x=temp.cg.flight,
                                                                   y=temp.cg.crack,
                                                                   xout=flt.last[temp.index])$y
                                
                                ksig.frst[temp.index] <- approx(x=temp.geo.crack,
                                                                y=temp.geo.ksig,
                                                                xout=a.frst[temp.index],
                                                                rule=2)$y
                                ksig.mddl[temp.index] <- approx(x=temp.geo.crack,
                                                                y=temp.geo.ksig,
                                                                xout=a.mddl[temp.index],
                                                                rule=2)$y
                                ksig.last[temp.index] <- approx(x=temp.geo.crack,
                                                                y=temp.geo.ksig,
                                                                xout=a.last[temp.index],
                                                                rule=2)$y
                            }

                        pos.frst <- pgumbel(kc, loc=(ms.lc*ksig.frst), scale=(ms.sc*ksig.frst))
                        pos.mddl <- pgumbel(kc, loc=(ms.lc*ksig.mddl), scale=(ms.sc*ksig.mddl))
                        pos.last <- pgumbel(kc, loc=(ms.lc*ksig.last), scale=(ms.sc*ksig.last))

                        ## simpson's rule approximation of the typical failure probability for
                        ##   each flight in the interval, taken to the power of the number
                        ##   of flights in the interval to estimate pr(fail) anywhere in the interval
                        pos.int   <- ( (pos.frst + 4*pos.mddl + pos.last) / 6 ) ^ num.flights

                        ## get weight update now and save it since this MAY POSSIBLY (not currently) be used to
                        ##   calculate sfpof.last assuming the only reason
                        ##   it could have failed during the interval is the K > Kc failure mode
                        ## w.kc.only <- w * pos.int
                        ## w.kc.only <- w.kc.only / sum(w.kc.only)

                        ##-------------------##
                        ## a>ac failure mode ##
                        ##-------------------##

                        ## particles which have breached the critical crack size
                        cc.hit             <- ( a.last >= cg.cc )
                        a.last[ cc.hit ]   <- Inf
                        flt.last[ cc.hit ] <- Inf
                        pos.int[ cc.hit ]  <- 0

                        ## probability of critical crack size breach during first or last flight:
                        pof.ac.frst <- rep(0, Np)
                        ## pof.ac.last <- pof.ac.frst
                        
                        ## simplest version: if failed anywhere in interval set equal probability to each flight
                        pof.ac.frst[cc.hit] <- 1 / num.flights
                        ## pof.ac.last[cc.hit] <- 1 / num.flights

                        pof.kc.frst <- 1 - pos.frst
                        ## pof.kc.last <- 1 - pos.last

                        ## not certain what this combination statement should be
                        pof.frst <- pof.ac.frst + pof.kc.frst - pof.ac.frst * pof.kc.frst
                        ## pof.last <- pof.ac.last + pof.kc.last - pof.ac.last * pof.kc.last

                        sfpof.first[iii] <- sum(w * pof.frst)
                        pof.int[iii]     <- sum(w * (1-pos.int))

                        ## bootstrap for SFPOF at first and last flights in the subinterval
                        ##   currently sfpof.last is just sfpof.frst, so use same boot for both also...
                        if( sfpof.boot.boolean )
                            {
                                sfpof.boot.n <- rep(0, boot.n)
                                Np.vec <- 1:Np
                                for(bbb in 1:boot.n)
                                    {
                                        index.boot <- sample(x=Np.vec, size=Np, replace=TRUE)
                                        pof.boot   <- pof.frst[index.boot]
                                        pof.boot   <- pof.frst[index.boot]
                                        w.boot     <- w[index.boot]
                                        w.boot     <- w.boot / sum(w.boot)

                                        sfpof.boot.n[bbb] <- sum(w.boot*pof.boot)

                                    }
                                sfpof.boot[2*iii-1,] <- quantile(x=sfpof.boot.n, probs=boot.q)
                                sfpof.boot[2*iii,]   <- sfpof.boot[2*iii-1,]
                            }

                        w <- w * pos.int
                        w <- w / sum(w)

                        ## XXX i would like to get a good solution for SFPOF of the last flight,
                        ##   but so far it is closer in general to use the estimate for the first flight only
                        ## sfpof.last[iii]  <- sum(w.kc.only * pof.last)
                        sfpof.last[iii]  <- sfpof.first[iii]
                        
                        a.frst   <- a.last
                        flt.frst <- flt.last

                        ## removed the silly flt.last + 1 stuff...see version 0.2-31 if you need it...

                    }

                obj$state$a <- a.frst
                obj$state$w <- w

            } else stop(paste(obj$parameters$survival.updating, "is not an option for parameter survival.updating"))

    sfpof.chronological <- as.numeric(rbind(sfpof.first, sfpof.last))

    if( sfpof.min > 0 )
        {
            sfpof.chronological[sfpof.chronological < sfpof.min] <- sfpof.min
            pof.int[pof.int < sfpof.min] <- sfpof.min
            if( sfpof.boot.boolean )
                {
                    sfpof.boot[sfpof.boot < sfpof.min] <- sfpof.min
                }
        }

    if( sfpof.boot.boolean )
        {
            obj$results$sfpof      <- rbind(obj$results$sfpof,
                                            data.frame(flight=calc.flights, sfpof=sfpof.chronological, sfpof.boot)
                                            )
        } else {
            obj$results$sfpof      <- rbind(obj$results$sfpof,
                                            data.frame(flight=calc.flights, sfpof=sfpof.chronological)
                                            )
        }

    obj$results$pof.int <- rbind(obj$results$pof.int,
                                 data.frame(flights.int = flights.per.calc,
                                            pof.int     = pof.int))
    
    return(obj)
}
