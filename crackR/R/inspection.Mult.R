inspection.Mult <-
function(obj, inspection.type=1)
{
    ## function for conducting a future inspection

    a   <- obj$state$a
    kc  <- obj$state$kc
    w   <- obj$state$w
    typ <- obj$state$typ

    ## if(any(is.na(w))) stop("all particles may have reached the critical crack length")
    if(any(is.na(w)))
        {
            ## temp debug
            stop("temporary error")
            warning("at least one weight is NaN. this usually means all particles have reached the critical crack length.")
            dta.types <- obj$parameters$dta.types
            repair.type.probs <- rep(NA, dta.types + 1)
            names(repair.type.probs) <- paste("pcd", 1:length(repair.type.probs), sep="")
            pcd.new <- data.frame(flight=max(obj$results$sfpof$flight),
                                  insp.type=inspection.type,
                                  t(repair.type.probs),
                                  pcd=NA)
            obj$results$pcd <- rbind(obj$results$pcd, pcd.new)
            return(obj)
        }

    bootstrap.pcd       <- obj$parameters$bootstrap.pcd
    bootstrap.samples   <- obj$parameters$bootstrap.samples
    bootstrap.quantiles <- obj$parameters$bootstrap.quantiles

    cg.cc         <- obj$parameters$cg.cc
    pod.threshold <- obj$parameters$pod.threshold
    dta.types     <- obj$parameters$dta.types

    dta <- obj$parameters$dta

    Np <- obj$parameters$Np

    pod.fun <- obj$parameters$pod.func[[inspection.type]]
    pod     <- pod.fun(a)

    pcd.vec <- pod*w
    pcd     <- sum(pcd.vec)

    ## need to loop through the repair types since each has it's own set of pod threshold values
    rep.types <- rep(0, Np)
    repair.type.probs <- rep(0, dta.types + 1)
    dta.types.current <- unique(typ)
    for(iii in dta.types.current)
        {
            thresholds.temp <- pod.threshold[iii,iii:dta.types] ## one line of the upper right triangle of threshold matrix
            index <- ( typ == iii )
            rep.types[index] <- findInterval(a[index], thresholds.temp) + iii ## yields repair type for each particle (if each were found)
        }
    rep.types.current <- unique(rep.types)
    for(jjj in rep.types.current) repair.type.probs[jjj] <- sum( pcd.vec[rep.types==jjj] )

    ## bootstrap PCD results
    if( bootstrap.pcd )
    {
        n.pcd <- dta.types + 2
        boot.results <- matrix(0, nrow=bootstrap.samples, ncol=n.pcd)
        
        for(bbb in 1:bootstrap.samples)
            {
                boot.index <- sample(1:Np, Np, replace=TRUE)
                boot.pcd.vec <- pcd.vec[ boot.index ]
                boot.rep.types <- rep.types[ boot.index ]

                boot.rep.type.probs <- rep(0, n.pcd - 1)
                boot.rep.types.current <- unique(boot.rep.types)
                for(jjj in boot.rep.types.current)
                    boot.rep.type.probs[jjj] <- sum( boot.pcd.vec[ boot.rep.types == jjj ] )
                boot.results[bbb,-n.pcd] <- boot.rep.type.probs
                boot.results[bbb,n.pcd]  <- sum( boot.pcd.vec )
                
            }
        boot.quantiles <- apply(boot.results, 2, function(x) quantile(x, bootstrap.quantiles))
    }

    ## downweight the particles according to the probability that they are found and repaired
    ## this will zero the weight of all particles at cg.cc (XXX why?! how?!)
    w <- w * (1-pod)
    w <- w / sum(w)  ## XXX should I be normalizing here ?!?! I don't think it's needed...

    ## remove particles at cg.cc since we don't want to re-sample these
    living     <- ( a < cg.cc )
    living.num <- sum(living)
    
    a   <- a[living]
    kc  <- kc[living]
    w   <- w[living]
    typ <- typ[living]
    
    ## set number of particles to come from new repairs and to keep from previous set
    ##   we want the total number of particles to get back to Np, with repaired count close to pcd*Np
    ##   also place a minimum number of particles in each repair type

    ## as of v02-14 the repair types include a last component of part replacement to type 1.
    ##   need to remove that last value and add it to the first value (type 1 repairs)
    pcd.by.type    <- repair.type.probs[-length(repair.type.probs)]
    pcd.by.type[1] <- pcd.by.type[1] + repair.type.probs[length(repair.type.probs)]

    ## first cut at number of particles to keep (Np.miss) and number for each repair type (Np.split)
    Np.miss  <- min( floor( Np * (1-pcd) ), living.num ) ## Np.miss can't be bigger than living.num
    Np.split <- floor( Np * pcd.by.type )

    ## increase repair particles to a minimum count
    Np.split.min <- 3
    Np.split[ Np.split < Np.split.min ] <- Np.split.min

    ## this block ensures there are a total of Np particles
    if( sum(Np.miss, Np.split) > Np )
        {
            ## if there are too many particles, drop the count for the largest category
            Np.largest <- which.max( c(Np.miss, Np.split) )
            if( Np.largest == 1 )
                {
                    Np.miss <- Np - sum( Np.split )
                } else {
                    Np.split[ (Np.largest-1) ] <- Np - Np.miss - sum( Np.split[ -(Np.largest-1) ] )
                }
        } else {
            ## if there are not enough particles, increase largest repair type
            ##   can't increase Np.miss because of living.num
            Np.split[ which.max( Np.split ) ] <- Np - Np.miss - sum( Np.split[ -which.max( Np.split ) ] )
        }

    ## missed (living) particles that are kept and normalized to weight (1-pcd)
    index <- sample(1:living.num, size=Np.miss, replace=FALSE)
    a  <- a[index]
    w  <- w[index]
    w  <- (1-pcd) * w / sum(w)
    kc <- kc[index]
    typ <- typ[index]

    ## gives integer values for the repair types that are being implemented in this set of repair particles
    ## if Np.split is "1" or a small number for a repair type, could this cause a problem
    ##   when the weights are normalized? do i need a minimum value for components of Np.split?
    ## dta.types.current <- (1:dta.types)[Np.split > 0]
    dta.types.current <- (1:dta.types) ## do this for every repair type!

    for(iii in dta.types.current)
        {
            a.rep   <- dta[[iii]]$ifs.rsamp(Np.split[iii])
            w.rep   <- dta[[iii]]$ifs.dactual(a.rep) / dta[[iii]]$ifs.dsamp(a.rep)

            ## new v0.2-35
            ##   just zero out the weights when this happens to avoid the possibility of an infinite loop
            ##   in this version, we count the number with weights greater than zero
            ##     if there are fewer than 3 (minimum particle count), we zero them all to prevent anomalies
            ## this can lead to sum(w) < 1, so needed to add the re-normalization before obj$state$a down below
            if( sum(w.rep > 0) < 3 )
                {
                    warning(paste("repair type", iii, "skipped, too few particles"))
                    w.rep <- rep(0, length(a.rep))
                } else {
                    w.rep   <- pcd.by.type[iii] * w.rep / sum(w.rep)
                }
            
            kc.rep  <- obj$parameters$dta[[iii]]$kc.rsamp(Np.split[iii])
            typ.rep <- rep(iii, Np.split[iii])
            
            a <- c(a, a.rep)
            w <- c(w, w.rep)
            kc <- c(kc, kc.rep)
            typ <- c(typ, typ.rep)
        }

    ## new v0.2-35 - see note above
    w <- w / sum(w)
    
    obj$state$a   <- a
    obj$state$w   <- w
    obj$state$kc  <- kc
    obj$state$typ <- typ

    names(repair.type.probs) <- paste("pcd", 1:length(repair.type.probs), sep="")

    pcd.new <- data.frame(flight=max(obj$results$sfpof$flight),
                          insp.type=inspection.type,
                          t(repair.type.probs),
                          pcd=pcd)  

    if( bootstrap.pcd )
    {
        pcd.boot <- as.data.frame(t(as.numeric(boot.quantiles)))
        temp.names <- c(paste("pcd", rep(1:(n.pcd-1), each=length(bootstrap.quantiles)), sep=""),
                        rep("pcd", each=length(bootstrap.quantiles)))
        names(pcd.boot) <- paste(temp.names, "_q",
                                 rep(bootstrap.quantiles, times=(n.pcd-1)), sep="")
        pcd.new <- data.frame(pcd.new, pcd.boot)
    }

    obj$results$pcd <- rbind(obj$results$pcd, pcd.new)

    return(obj)
}
