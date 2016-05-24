# from parallel:::clusterSetRNGStream(GPL-2)
Rhpc_setupRNG <- function(cl, iseed = NULL)
{
    oldseed <-
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        else NULL
    RNGkind("L'Ecuyer-CMRG")
    if(!is.null(iseed)) set.seed(iseed)
    nc <- Rhpc_numberOfWorker(cl)
    seeds <- vector("list", nc)
    seeds[[1L]] <- .Random.seed
    for(i in seq_len(nc-1L)) seeds[[i+1L]] <- parallel::nextRNGStream(seeds[[i]])
    ## Reset the random seed in the master.
    if(!is.null(oldseed))
        assign(".Random.seed", oldseed, envir = .GlobalEnv)
    else rm(.Random.seed, envir = .GlobalEnv)
    Rhpc_lapply(cl, seeds, function(seed) assign(".Random.seed", seed, envir = .GlobalEnv))
    invisible()
}


