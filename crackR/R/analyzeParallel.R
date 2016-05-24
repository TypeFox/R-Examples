analyzeParallel <-
function(obj, n.parallel = 3)
{
    ## simply runs analyze sequentially. intended that later it could use parallel processing
    ## no "add" option exists at present; not sure how that would make sense
    ##   could keep the same number of sequences and just beef each up with samples

    ## if n.parallel < 2, stop
    if( n.parallel < 2 ) stop("n.parallel should be at least 2")
    
    ## there's no "add" option, so if a crackRparallel run is submitted just use first parameter set
    if( any(class(obj) == "crackRparallel") ) obj <- obj[[1]]$parameters
    
    ## since running several times, if it's already initialized just drop back to the parameter set
    if( any(class(obj) == "Sing") || any(class(obj) == "Mult") || any(class(obj) == "CD") ) obj <- obj$parameters
    
    ## run analyze repeatedly, storing each result
    parallel.list <- list()
    length(parallel.list) <- n.parallel
    for(iii in 1:n.parallel)
        {
            cat("Beginning parallel run #", iii, "\n")
            temp.init <- crackRinit(obj)
            print(system.time(parallel.list[[iii]] <- analyze(temp.init, add=FALSE)$results))
        }

    ## adding a new component that houses the combined results
    combined.sfpof   <- parallel.list[[1]]$sfpof$sfpof
    combined.pof.int <- parallel.list[[1]]$pof.int$pof.int
    combined.pcd     <- parallel.list[[1]]$pcd[,-c(1:2)]
    for(iii in 2:n.parallel)
    {
        combined.sfpof   <- combined.sfpof   + parallel.list[[iii]]$sfpof$sfpof
        combined.pof.int <- combined.pof.int + parallel.list[[iii]]$pof.int$pof.int
        combined.pcd     <- combined.pcd     + parallel.list[[iii]]$pcd[,-c(1:2)]
    }
    combined <- parallel.list[[1]]
    combined$Np.total <- n.parallel * parallel.list[[1]]$Np.total
    combined$sfpof$sfpof     <- combined.sfpof   / n.parallel
    combined$pof.int$pof.int <- combined.pof.int / n.parallel
    combined$pcd[,-c(1:2)]   <- combined.pcd     / n.parallel

    parallel.out <- list(parallel.list, combined)
    class(parallel.out) <-  c("crackRparallel", "list")

    return(parallel.out)
}
