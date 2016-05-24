calcInterval.CD <-
function(obj, interval.flights)
{
    ## takes a crackR CD object and advances a specified interval of flights
    
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

    dta.pc <- obj$parameters$dta.pc ## deterministic dta data - primary cold
    dta.ph <- obj$parameters$dta.ph ## deterministic dta data - primary hot
    dta.sc <- obj$parameters$dta.sc ## deterministic dta data - secondary cold
    dta.sh <- obj$parameters$dta.sh ## deterministic dta data - secondary hot
    
    ms.lc     <- obj$parameters$ms.gumbel[1] ## gumbel max stress dist is hard coded at present
    ms.sc     <- obj$parameters$ms.gumbel[2] ## " "
    sfpof.min <- obj$parameters$sfpof.min    ## smaller estimates of sfpof are truncated at this value

    cg.cc.pc <- obj$parameters$cg.cc.pc ## critical crack length for primary at which secondary goes hot
    cg.cc.ph <- obj$parameters$cg.cc.ph ## critical crack length at which pr(fail)=100%, primary crack
    cg.cc.sc <- obj$parameters$cg.cc.sc ## critical crack length for secondary at which primary goes hot
    cg.cc.sh <- obj$parameters$cg.cc.sh ## critical crack length at which pr(fail)=100%, secondary

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

    ## store the proportion of particles that are hot for cont. dam. diagnostics
    p.cold.frst  <- rep(0, calc.num)
    p.cold.last  <- rep(0, calc.num)
    s.cold.frst  <- rep(0, calc.num)
    s.cold.last  <- rep(0, calc.num)

    ## model state at the start of the routine
    a.p <- obj$state$a.p
    a.s <- obj$state$a.s
    kc  <- obj$state$kc
    w   <- obj$state$w

    Np  <- obj$parameters$Np

    ## split calculation here for flight-by-flight or interval calculation

    if( tolower(obj$parameters$survival.updating) == "fbf" )
        {

            ## initialize current flight numbers and ksig values for later masking
            flt.p  <- rep(0, Np)
            flt.s  <- flt.p
            ksig.p <- flt.p
            ksig.s <- flt.p

            ## identify initial cold/hot status of each particle (will be repeatedly overwritten)
            primary.cold   <- ( a.s < cg.cc.sc ) ## secondary is not yet failed
            secondary.cold <- ( a.p < cg.cc.pc ) ##   primary is not yet failed

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
                                stop("at least one weight is NaN. this usually means all particles have reached the critical crack length.")
                            }

                    sfpof.temp <- rep(0, flights.per.calc[iii])
                    
                    p.cold.frst[iii] <- sum(   primary.cold * w)
                    s.cold.frst[iii] <- sum( secondary.cold * w)

                    for(jjj in 1:flights.per.calc[iii])
                        {

                            ## re-reference the hot/cold status
                            primary.cold   <- ( a.s < cg.cc.sc )
                            secondary.cold <- ( a.p < cg.cc.pc )

                            ## indices for growing cracks and finding stress intensities
                            ##   XXX might be more efficient to remove cracks of length Inf (likely not much)
                            index.pc <- primary.cold
                            index.ph <- !primary.cold
                            index.sc <- secondary.cold
                            index.sh <- !secondary.cold

                            ## reference the first flight using the crack sizes and the hot/cold status
                            ## this must be done in each calculation loop because the status
                            ##   may switch from cold to hot (different current flight numbers)
                            flt.p[index.pc] <- approxExtrap(x=dta.pc$cg$crack,
                                                            y=dta.pc$cg$flight,
                                                            xout=a.p[index.pc])$y
                            flt.p[index.ph] <- approxExtrap(x=dta.ph$cg$crack,
                                                            y=dta.ph$cg$flight,
                                                            xout=a.p[index.ph])$y
                            flt.s[index.sc] <- approxExtrap(x=dta.sc$cg$crack,
                                                            y=dta.sc$cg$flight,
                                                            xout=a.s[index.sc])$y
                            flt.s[index.sh] <- approxExtrap(x=dta.sh$cg$crack,
                                                            y=dta.sh$cg$flight,
                                                            xout=a.s[index.sh])$y

                            ## crack sizes this flight
                            a.p[index.pc] <- approxExtrap(x=dta.pc$cg$flight,
                                                          y=dta.pc$cg$crack,
                                                          xout=flt.p[index.pc]+1)$y
                            a.p[index.ph] <- approxExtrap(x=dta.ph$cg$flight,
                                                          y=dta.ph$cg$crack,
                                                          xout=flt.p[index.ph]+1)$y
                            a.s[index.sc] <- approxExtrap(x=dta.sc$cg$flight,
                                                          y=dta.sc$cg$crack,
                                                          xout=flt.s[index.sc]+1)$y
                            a.s[index.sh] <- approxExtrap(x=dta.sh$cg$flight,
                                                          y=dta.sh$cg$crack,
                                                          xout=flt.s[index.sh]+1)$y

                            ## K/sigma values this flight
                            ksig.p[index.pc] <- approx(x=dta.pc$geo$crack,
                                                       y=dta.pc$geo$ksig,
                                                       xout=a.p[index.pc],
                                                       rule=2)$y
                            ksig.p[index.ph] <- approx(x=dta.ph$geo$crack,
                                                       y=dta.ph$geo$ksig,
                                                       xout=a.p[index.ph],
                                                       rule=2)$y
                            ksig.s[index.sc] <- approx(x=dta.sc$geo$crack,
                                                       y=dta.sc$geo$ksig,
                                                       xout=a.s[index.sc],
                                                       rule=2)$y
                            ksig.s[index.sh] <- approx(x=dta.sh$geo$crack,
                                                       y=dta.sh$geo$ksig,
                                                       xout=a.s[index.sh],
                                                       rule=2)$y
                            
                            ## failure probability
                            ##   for continuing damage the crack with LARGER survival probability
                            ##     drives the SFPOF estimates since both must fail (due to same load)

                            pos.current.p  <- pgumbel(kc, loc=(ms.lc*ksig.p), scale=(ms.sc*ksig.p)) 
                            pos.current.s  <- pgumbel(kc, loc=(ms.lc*ksig.s), scale=(ms.sc*ksig.s)) 

                            ## a>ac failure mode
                            cc.hit.p <- ( a.p >= cg.cc.ph )
                            cc.hit.s <- ( a.s >= cg.cc.sh )
                            pos.current.p[ cc.hit.p ] <- 0 ## a>ac failure mode
                            pos.current.s[ cc.hit.s ] <- 0

                            ## cracks larger than critical size set to Inf so easily identified
                            a.p[ cc.hit.p ] <- Inf
                            a.s[ cc.hit.s ] <- Inf

                            ## for cont dam, the relevant survival probability is the larger one (same applied load)
                            pos.current <- apply(cbind(pos.current.p, pos.current.s), 1, max)
                            
                            sfpof.temp[jjj] <- 1 - sum(w * pos.current)

                            ## bootstrap for SFPOF at first and last flight in the subinterval
                            if( sfpof.boot.boolean )
                            {
                                if( ( jjj == 1 ) || ( jjj == flights.per.calc[iii] ) )
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
                            }

                            w <- w * pos.current
                            w <- w / sum(w)   
                        }
                    
                    sfpof.first[iii] <- sfpof.temp[1]
                    sfpof.last[iii]  <- sfpof.temp[flights.per.calc[iii]]
                    pof.int[iii]     <- 1-prod(1-sfpof.temp)

                    p.cold.last[iii] <- sum(   primary.cold * w)
                    s.cold.last[iii] <- sum( secondary.cold * w)

                }

            obj$state$a.p <- a.p
            obj$state$a.s <- a.s
            obj$state$w <- w

        } else if( tolower(obj$parameters$survival.updating) == "int" )
            {

                a.p.frst <- a.p
                a.s.frst <- a.s

                ## initialize the various vectors that are filled with masking in loop
                flt.p.frst      <- rep(0, Np)
                flt.p.mddl      <- flt.p.frst
                flt.p.last      <- flt.p.frst
                a.p.mddl        <- flt.p.frst
                a.p.last        <- flt.p.frst
                ksig.p.frst     <- flt.p.frst
                ksig.p.mddl     <- flt.p.frst
                ksig.p.last     <- flt.p.frst
                flt.s.frst      <- flt.p.frst
                flt.s.mddl      <- flt.p.frst
                flt.s.last      <- flt.p.frst
                a.s.mddl        <- flt.p.frst
                a.s.last        <- flt.p.frst
                ksig.s.frst     <- flt.p.frst
                ksig.s.mddl     <- flt.p.frst
                ksig.s.last     <- flt.p.frst
                
                ## identify the cold/hot status of each particle
                primary.cold   <- ( a.s.frst < cg.cc.sc ) ## secondary is not yet failed
                secondary.cold <- ( a.p.frst < cg.cc.pc ) ##   primary is not yet failed

                ## set flt.p and flt.s according to type and cold/hot status AT THE START OF THE INTERVAL
                index.pc <- primary.cold
                index.ph <- !primary.cold
                index.sc <- secondary.cold
                index.sh <- !secondary.cold

                ## using approx so we don't start with hugely negative flight hours
                flt.p.frst[index.pc] <- approx(x=dta.pc$cg$crack,
                                               y=dta.pc$cg$flight,
                                               xout=a.p.frst[index.pc],
                                               rule=2)$y
                flt.p.frst[index.ph] <- approx(x=dta.ph$cg$crack,
                                               y=dta.ph$cg$flight,
                                               xout=a.p.frst[index.ph],
                                               rule=2)$y
                flt.s.frst[index.sc] <- approx(x=dta.sc$cg$crack,
                                               y=dta.sc$cg$flight,
                                               xout=a.s.frst[index.sc],
                                               rule=2)$y
                flt.s.frst[index.sh] <- approx(x=dta.sh$cg$crack,
                                               y=dta.sh$cg$flight,
                                               xout=a.s.frst[index.sh],
                                               rule=2)$y

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
                                    stop("at least one weight is NaN. this usually means all particles have reached the critical crack length.")
                                }
                        
                        ## frst = first flight in the interval
                        ## mddl = middle flight in int
                        ## last = last flight in int

                        num.flights <- flights.per.calc[iii]

                        ##-------------------##
                        ## K>Kc failure mode ##
                        ##-------------------##

                        ## hot/cold status
                        ## here could check if flt.p/s.last is off the table...
                        primary.cold   <- ( a.s.frst < cg.cc.sc ) ## secondary is not yet failed
                        secondary.cold <- ( a.p.frst < cg.cc.pc ) ##   primary is not yet failed

                        p.cold.frst[iii] <- sum(   primary.cold * w)
                        s.cold.frst[iii] <- sum( secondary.cold * w)

                        ## might be faster to un-index any cracks of size Inf
                        index.pc <- primary.cold
                        index.ph <- !primary.cold
                        index.sc <- secondary.cold
                        index.sh <- !secondary.cold

                        ## flt.p/s.frst needs re-referenced in case cold/hot status changed
                        flt.p.frst[index.pc] <- approxExtrap(x=dta.pc$cg$crack,
                                                             y=dta.pc$cg$flight,
                                                             xout=a.p.frst[index.pc])$y
                        flt.p.frst[index.ph] <- approxExtrap(x=dta.ph$cg$crack,
                                                             y=dta.ph$cg$flight,
                                                             xout=a.p.frst[index.ph])$y
                        flt.s.frst[index.sc] <- approxExtrap(x=dta.sc$cg$crack,
                                                             y=dta.sc$cg$flight,
                                                             xout=a.s.frst[index.sc])$y
                        flt.s.frst[index.sh] <- approxExtrap(x=dta.sh$cg$crack,
                                                             y=dta.sh$cg$flight,
                                                             xout=a.s.frst[index.sh])$y

                        flt.p.mddl <- flt.p.frst + num.flights*0.5 + 0.5 ## mddl of mddl flight
                        flt.s.mddl <- flt.s.frst + num.flights*0.5 + 0.5
                        flt.p.last <- flt.p.frst + num.flights + 1 ## end of last flight
                        flt.s.last <- flt.s.frst + num.flights + 1

                        a.p.mddl[index.pc] <- approxExtrap(x=dta.pc$cg$flight,
                                                           y=dta.pc$cg$crack,
                                                           xout=flt.p.mddl[index.pc])$y
                        a.p.mddl[index.ph] <- approxExtrap(x=dta.ph$cg$flight,
                                                           y=dta.ph$cg$crack,
                                                           xout=flt.p.mddl[index.ph])$y
                        a.s.mddl[index.sc] <- approxExtrap(x=dta.sc$cg$flight,
                                                           y=dta.sc$cg$crack,
                                                           xout=flt.s.mddl[index.sc])$y
                        a.s.mddl[index.sh] <- approxExtrap(x=dta.sh$cg$flight,
                                                           y=dta.sh$cg$crack,
                                                           xout=flt.s.mddl[index.sh])$y

                        a.p.last[index.pc] <- approxExtrap(x=dta.pc$cg$flight,
                                                           y=dta.pc$cg$crack,
                                                           xout=flt.p.last[index.pc])$y
                        a.p.last[index.ph] <- approxExtrap(x=dta.ph$cg$flight,
                                                           y=dta.ph$cg$crack,
                                                           xout=flt.p.last[index.ph])$y
                        a.s.last[index.sc] <- approxExtrap(x=dta.sc$cg$flight,
                                                           y=dta.sc$cg$crack,
                                                           xout=flt.s.last[index.sc])$y
                        a.s.last[index.sh] <- approxExtrap(x=dta.sh$cg$flight,
                                                           y=dta.sh$cg$crack,
                                                           xout=flt.s.last[index.sh])$y

                        ksig.p.frst[index.pc] <- approx(x=dta.pc$geo$crack,
                                                        y=dta.pc$geo$ksig,
                                                        xout=a.p.frst[index.pc], rule=2)$y
                        ksig.p.frst[index.ph] <- approx(x=dta.ph$geo$crack,
                                                        y=dta.ph$geo$ksig,
                                                        xout=a.p.frst[index.ph], rule=2)$y
                        ksig.s.frst[index.sc] <- approx(x=dta.sc$geo$crack,
                                                        y=dta.sc$geo$ksig,
                                                        xout=a.s.frst[index.sc], rule=2)$y
                        ksig.s.frst[index.sh] <- approx(x=dta.sh$geo$crack,
                                                        y=dta.sh$geo$ksig,
                                                        xout=a.s.frst[index.sh], rule=2)$y

                        ksig.p.mddl[index.pc] <- approx(x=dta.pc$geo$crack,
                                                        y=dta.pc$geo$ksig,
                                                        xout=a.p.mddl[index.pc], rule=2)$y
                        ksig.p.mddl[index.ph] <- approx(x=dta.ph$geo$crack,
                                                        y=dta.ph$geo$ksig,
                                                        xout=a.p.mddl[index.ph], rule=2)$y
                        ksig.s.mddl[index.sc] <- approx(x=dta.sc$geo$crack,
                                                        y=dta.sc$geo$ksig,
                                                        xout=a.s.mddl[index.sc], rule=2)$y
                        ksig.s.mddl[index.sh] <- approx(x=dta.sh$geo$crack,
                                                        y=dta.sh$geo$ksig,
                                                        xout=a.s.mddl[index.sh], rule=2)$y

                        ksig.p.last[index.pc] <- approx(x=dta.pc$geo$crack,
                                                        y=dta.pc$geo$ksig,
                                                        xout=a.p.last[index.pc], rule=2)$y
                        ksig.p.last[index.ph] <- approx(x=dta.ph$geo$crack,
                                                        y=dta.ph$geo$ksig,
                                                        xout=a.p.last[index.ph], rule=2)$y
                        ksig.s.last[index.sc] <- approx(x=dta.sc$geo$crack,
                                                        y=dta.sc$geo$ksig,
                                                        xout=a.s.last[index.sc], rule=2)$y
                        ksig.s.last[index.sh] <- approx(x=dta.sh$geo$crack,
                                                        y=dta.sh$geo$ksig,
                                                        xout=a.s.last[index.sh], rule=2)$y

                        pos.p.frst <- pgumbel(kc, loc=(ms.lc*ksig.p.frst), scale=(ms.sc*ksig.p.frst))
                        pos.p.mddl <- pgumbel(kc, loc=(ms.lc*ksig.p.mddl), scale=(ms.sc*ksig.p.mddl))
                        pos.p.last <- pgumbel(kc, loc=(ms.lc*ksig.p.last), scale=(ms.sc*ksig.p.last))

                        pos.s.frst <- pgumbel(kc, loc=(ms.lc*ksig.s.frst), scale=(ms.sc*ksig.s.frst))
                        pos.s.mddl <- pgumbel(kc, loc=(ms.lc*ksig.s.mddl), scale=(ms.sc*ksig.s.mddl))
                        pos.s.last <- pgumbel(kc, loc=(ms.lc*ksig.s.last), scale=(ms.sc*ksig.s.last))

                        ## simpson's rule
                        pos.p.int   <- ( (pos.p.frst + 4*pos.p.mddl + pos.p.last) / 6 ) ^ num.flights
                        pos.s.int   <- ( (pos.s.frst + 4*pos.s.mddl + pos.s.last) / 6 ) ^ num.flights

                        #######################
                        ## a>ac failure mode ##
                        #######################

                        cc.hit.p <- ( a.p.last >= cg.cc.ph )
                        cc.hit.s <- ( a.s.last >= cg.cc.sh )
                        pos.p.int[cc.hit.p] <- 0
                        pos.s.int[cc.hit.s] <- 0

                        ## for cont dam the LARGER survival probability drives SFPOF updating
                        pos.int <- apply(cbind(pos.p.int, pos.s.int), 1, max)

                        ## as in the ordinary fast version, sfpof is estimated a bit differently here
                        ## we use ALL the critical crack failures to help the estimate at the
                        ##   first and last flights in the interval

                        ## probability of critical crack size breach during first or last flight:
                        pof.ac.p <- rep(0, Np)
                        pof.ac.s <- pof.ac.p
                        pof.ac.p[cc.hit.p] <- 1 / num.flights
                        pof.ac.s[cc.hit.s] <- 1 / num.flights

                        pof.kc.p.frst <- 1 - pos.p.frst
                        pof.kc.s.frst <- 1 - pos.s.frst

                        ## not certain what this combination statement should be
                        pof.p.frst <- pof.ac.p + pof.kc.p.frst - pof.ac.p * pof.kc.p.frst
                        pof.s.frst <- pof.ac.s + pof.kc.s.frst - pof.ac.s * pof.kc.s.frst

                        ## cont dam, use LOWER failure probability for sfpof estimate
                        pof.frst <- apply(cbind(pof.p.frst, pof.s.frst), 1, min)
                        
                        sfpof.first[iii] <- sum(w * pof.frst)
                        pof.int[iii]     <- sum(w * (1-pos.int))

                        ## cg.cc cracks set to size Inf for easy identification
                        a.p.last[cc.hit.p] <- Inf
                        a.s.last[cc.hit.s] <- Inf
                        
                        a.p.frst <- a.p.last
                        a.s.frst <- a.s.last

                        p.cold.last[iii] <- sum( ( a.s.last < cg.cc.sc ) * w)
                        s.cold.last[iii] <- sum( ( a.p.last < cg.cc.pc ) * w)

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

                        ## XXX i would like to get a good solution for the SFPOF of the last flight,
                        ##   but so far it is closer in general to use the estimate for the first flight only
                        sfpof.last[iii]  <- sfpof.first[iii]
                        
                    }

                obj$state$a.p <- a.p.last
                obj$state$a.s <- a.s.last
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

    pri.cold.prop <- as.numeric(rbind(p.cold.frst, p.cold.last))
    sec.cold.prop <- as.numeric(rbind(s.cold.frst, s.cold.last))
    
    if( sfpof.boot.boolean )
    {
        obj$results$sfpof      <- rbind(obj$results$sfpof,
                                        data.frame(flight=calc.flights, sfpof=sfpof.chronological,
                                                   pri.cold.prop=pri.cold.prop, sec.cold.prop=sec.cold.prop,
                                                   sfpof.boot=sfpof.boot)
                                        )
    } else {
        obj$results$sfpof      <- rbind(obj$results$sfpof,
                                        data.frame(flight=calc.flights, sfpof=sfpof.chronological,
                                                   pri.cold.prop=pri.cold.prop, sec.cold.prop=sec.cold.prop)
                                        )
    }

    obj$results$pof.int <- rbind(obj$results$pof.int,
                                 data.frame(flights.int = flights.per.calc,
                                            pof.int     = pof.int))

    return(obj)
}
