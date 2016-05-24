# parallel version of lapply, using parallel::mclapply. We loose verbosity.

parallel_lapply <- function(X,FUN,...,mc.cores,mc.set.seed=FALSE,mc.silent=TRUE,verbose=TRUE)
{
    if (!is.vector(X) || is.object(X)) 
        X <- as.list(X)
    last_disp <- 0
    cat("\r")
    sx <- seq_along(X)
    res <- vector("list", length(sx))
        names(res) <- names(X)

    cores <- as.integer(mc.cores)

    if(.Platform$OS.type=="unix" && cores > 1)
    {
        do_parallel <- TRUE
    }
    else
    {
        do_parallel <- FALSE
    }

    if(do_parallel)
    {
        if(verbose)
        {
            cat("\r")
            message <- paste(
                'Executing ',
                length(sx),
                ' jobs in parallel',
                    sep=''
                )
            last_disp <- nchar(message)
            cat(message)
        }
        res <- parallel::mclapply(X,FUN,...,mc.cores=mc.cores,mc.set.seed=FALSE,mc.silent=TRUE)
    }
    else
    {
        k<-0
        nb_jobs_total <- length(sx)
        nb_jobs_done <- 0
        for(x in sx)
        {
            if(verbose)
            {
                cat("\r")
                for(theta in 1:last_disp)
                {
                    cat(" ")
                }
                cat("\r")
                message <- paste(
                    'Non-parallel jobs: ',
                    nb_jobs_done,
                    '/',
                    nb_jobs_total,
                    ' ',
                    sep=''
                )
                last_disp <- nchar(message)
                cat(message)
            }
            
            res[[x]] <- FUN(X[[x]],...)
            nb_jobs_done <- nb_jobs_done + 1
        }
    }
    
    if(verbose)
    {
        cat("\r")    
        for(theta in 1:last_disp)
        {
            cat(" ")
        }
        cat("\r")
    }


    return(res)
}



