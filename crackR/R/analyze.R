analyze <-
function(obj, add=FALSE)
{
    ## runs a crackR analysis using a previously initialized crackR object...
    ## can be run using a new crackR object or one that contains previous results
    ## to add to the existing results, select add=TRUE
    
    ## if parameters are submitted, try initializing it
    if( class(obj)[1] == "list" )
    {
        obj <- crackRinit(obj)
        warning("uninitialized object submitted; attempting to initialize.")
    }
    
    ## if adding results to an existing run, save previous information before running a fresh analysis
    if(add==TRUE)
        {
            previous.results <- obj$results
            new.Np.total     <- obj$parameters$Np + previous.results$Np.total

            obj <- crackRinit(obj$parameters)
            
        } else {
            ## if we are not adding to existing results, see if there are any..if so, overwrite
            if(length(obj$results$sfpof$sfpof) > 0) obj <- crackRinit(obj$parameters)
        }

    ## run analysis one interval at a time
    insp <- obj$parameters$inspections
    for(jjj in 1:dim(insp)[1])
        {
            obj <- calcInterval(obj, insp$flt.interval[jjj])
            obj <- inspection(obj, inspection.type=insp$type[jjj])
        }

    ## if adding results to an existing run, combine new results with previous results as additional results component
    ## save new results in a "new.results" list
    if(add==TRUE)
        {
            obj$new.results <- obj$results
            new.run.prop    <- obj$parameters$Np / new.Np.total

            flight <- obj$results$sfpof$flight

            sfpof.df  <- new.run.prop * obj$results$sfpof + (1-new.run.prop) * previous.results$sfpof ## also factoring flight!
            sfpof.df$flight <- flight ## doing this to preserve the data frame structure of results$sfpof

            flights.int <- obj$results$pof.int$flights.int
            pof.int <- new.run.prop * obj$results$pof.int$pof.int + (1-new.run.prop) * previous.results$pof.int$pof.int
            pof.int <- data.frame(flights.int=flights.int, pof.int=pof.int)

            insp.flight <- obj$results$pcd$flight
            insp.type   <- obj$parameters$inspections$type

            pcd.df     <- obj$results$pcd[,-(1:2)]
            pcd.new.df <- previous.results$pcd[,-(1:2)]
            pcd.out    <- new.run.prop * pcd.df  + (1-new.run.prop) * pcd.new.df

            pcd  <- data.frame(flight=insp.flight, insp.type=insp.type, pcd.out)
            
            obj$results     <- list(Np.total = new.Np.total, sfpof = sfpof.df, pof.int=pof.int, pcd = pcd)

            class(obj$results) <- c("crackRresults", "list")
        }
    
    return(obj)
}
