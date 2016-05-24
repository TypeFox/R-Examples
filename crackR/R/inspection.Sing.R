inspection.Sing <-
function(obj, inspection.type=1)
{
    ## function for conducting a future inspection with a single DTA type
    ## as of version 0.3-03, bootstrapping PCD is an option

    a   <- obj$state$a
    kc  <- obj$state$kc
    w   <- obj$state$w

    ## if(any(is.na(w))) stop("all particles may have reached the critical crack length")
    if(any(is.na(w)))
        stop("at least one weight is NaN. this usually means all particles have reached the critical crack length.")

    bootstrap.pcd       <- obj$parameters$bootstrap.pcd
    bootstrap.samples   <- obj$parameters$bootstrap.samples
    bootstrap.quantiles <- obj$parameters$bootstrap.quantiles
    
    cg.cc <- obj$parameters$cg.cc

    pod.threshold <- obj$parameters$pod.threshold
    ## for single DTA type analysis this should be a vector
    if( class( pod.threshold ) == "matrix" ) stop("pod.threshold should be a vector for analysis type \"single\"")

    dta <- obj$parameters$dta

    Np <- obj$parameters$Np

    pod.fun <- obj$parameters$pod.func[[inspection.type]]
    pod     <- pod.fun(a)

    pcd.vec <- pod*w
    pcd     <- sum(pcd.vec)

    repair.type.probs <- rep(0, length(pod.threshold) + 1)

    ## yields repair type for each particle (if each were found)
    rep.types <- findInterval(a, pod.threshold) + 1

    ## what repair types exist here? (for looping)
    rep.types.current <- unique(rep.types)
    for(jjj in rep.types.current)
        repair.type.probs[jjj] <- sum( pcd.vec[ rep.types == jjj ] )
    ## version 0.3-04, removed the first component of the sum since ALWAYS zero for Sing case

    ## bootstrap PCD results
    if( bootstrap.pcd )
    {
        n.pcd <- length(pod.threshold) + 2
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
    w <- w * (1-pod)
    w <- w / sum(w)

    ## remove particles at cg.cc since we don't want to re-sample these
    living     <- ( a < cg.cc )
    living.num <- sum(living)
    
    a   <- a[living]
    kc  <- kc[living]
    w   <- w[living]
    
    ## set number of particles to come from new repairs and to keep from previous set
    ##   we want the total number of particles to get back to Np, with repaired count close to pcd*Np

    ## the repair types include a last component of part replacement to type 1.
    ##   need to remove that last value and add it to the first value (type 1 repairs)
    pcd.by.type    <- repair.type.probs[-length(repair.type.probs)]
    pcd.by.type[1] <- pcd.by.type[1] + repair.type.probs[length(repair.type.probs)]

    ## first cut at number of particles to keep (Np.miss) and number to repair (Np.rep)
    Np.miss <- min( floor( Np * (1-pcd) ), living.num ) ## Np.miss can't be bigger than living.num
    Np.rep  <- Np - Np.miss

    ## increase particle counts to a minimum value
    Np.min <- 10
    if( Np.miss < Np.min ) Np.miss <- Np.min
    if( Np.rep  < Np.min ) Np.rep  <- Np.min
    
    ## this block ensures there are a total of Np particles
    if( sum(Np.miss, Np.rep) > Np )
        {
            ## there are too many particles, so drop the larger count
            if( Np.miss > Np.rep )
                {
                    Np.miss <- Np - Np.rep
                } else {
                    Np.rep  <- Np - Np.miss
                }
        } else {
            ## if there are not enough particles, increase repair count
            ##   can't increase Np.miss because of living.num
            Np.rep <- Np - Np.miss
        }

    ## missed (living) particles that are kept and normalized to weight (1-pcd)
    index <- sample(1:living.num, size=Np.miss, replace=FALSE)
    a  <- a[index]
    w  <- w[index]
    w  <- (1-pcd) * w / sum(w)
    kc <- kc[index]

    a.rep <- dta$rfs.rsamp(Np.rep)
    w.rep <- dta$rfs.dactual(a.rep) / dta$rfs.dsamp(a.rep)

    ## if all (or almost all) repair particles have zero weight, this wreaks havoc
    ## make sure at least 3 have positive weight, otherwise zero them all
    if( sum(w.rep > 0) < 3 )
        {
            warning(paste("repair skipped, too few particles"))
            w.rep <- rep(0, length(a.rep))
        } else {
            w.rep   <- pcd * w.rep / sum(w.rep)
        }
            
    kc.rep  <- obj$parameters$dta$kc.rsamp(Np.rep)
    
    a <- c(a, a.rep)
    w <- c(w, w.rep)
    kc <- c(kc, kc.rep)

    ## re-normalize the weights in case the repair weights were zeroed out above
    w <- w / sum(w)
    
    obj$state$a   <- a
    obj$state$w   <- w
    obj$state$kc  <- kc

    names(repair.type.probs) <- paste("pcd", 1:length(repair.type.probs), sep="")

    pcd.new <- data.frame(flight=max(obj$results$sfpof$flight),
                          insp.type=inspection.type,
                          t(repair.type.probs),
                          pcd=pcd)  

    if( bootstrap.pcd )
    {
        pcd.boot <- as.data.frame(t(as.numeric(boot.quantiles)))
        temp.names <- c(paste("pcd", rep(1:(length(pod.threshold)+1), each=length(bootstrap.quantiles)), sep=""),
                        rep("pcd", each=length(bootstrap.quantiles)))
        names(pcd.boot) <- paste(temp.names, "_q",
                                rep(bootstrap.quantiles, times=(length(pod.threshold)+1)), sep="")
        pcd.new <- data.frame(pcd.new, pcd.boot)
    }
    
    obj$results$pcd <- rbind(obj$results$pcd, pcd.new)

    return(obj)
}
