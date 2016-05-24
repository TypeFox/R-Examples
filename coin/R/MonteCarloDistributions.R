split_index <- function(n, by) {
    if (n < by)
        by <- n
    lengths(lapply(seq_len(by), function(i) seq.int(i, n, by)),
            use.names = FALSE)
}

MonteCarlo <- function(x, y, block, weights, B, parallel, ncpus, cl) {
    ## expand observations for non-unit weights
    if (!is_unity(weights)) {
        idx <- rep.int(seq_along(weights), weights)
        x <- x[idx, , drop = FALSE]
        y <- y[idx, , drop = FALSE]
        block <- block[idx]
    }

    montecarlo <- function(B)
        .Call("R_MonteCarloIndependenceTest",
              x, y, as.integer(block), as.integer(B), PACKAGE = "coin")

    if (parallel == "no")
        montecarlo(B)
    else {
        ## load the 'parallel' namespace if necessary
        if (!isNamespaceLoaded("parallel")) {
            ## loading 'parallel' changes RNG state if R_PARALLEL_PORT is unset
            if (Sys.getenv("R_PARALLEL_PORT") == "") {
                ## make sure '.Random.seed' exists; almost unnecessary since it
                ## always does when called from "ApproxNullDistribution"
                if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
                    runif(1L)
                ## save existing RNG state
                seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
                ## load namespace
                if (!requireNamespace("parallel", quietly = TRUE))
                    stop("package ", sQuote("parallel"),
                         " is needed for parallel operation")
                ## put back the saved RNG state
                assign(".Random.seed", value = seed, envir = .GlobalEnv)
            } else
                ## load namespace
                if (!requireNamespace("parallel", quietly = TRUE))
                    stop("package ", sQuote("parallel"),
                         " is needed for parallel operation")
        }

        if (RNGkind()[1L] == "L'Ecuyer-CMRG")
            ## advance stream in master process upon exit
            on.exit(assign(".Random.seed",
                           value = parallel::nextRNGStream(
                               get(".Random.seed", envir = .GlobalEnv,
                                   inherits = FALSE)),
                           envir = .GlobalEnv))

        if (parallel == "multicore") {
            if (.Platform$OS.type == "windows")
                stop(sQuote(paste0("parallel = ", dQuote("multicore"))),
                     " is not available for MS Windows")
            if (as.integer(ncpus) < 2L)
                warning("parallel operation requires at least two processes")
###            Bp <- split_index(B, ncpus) # distribute workload evenly
###            RET <- parallel::mclapply(Bp, FUN = montecarlo, mc.cores = ncpus,
###                                      mc.allow.recursive = FALSE)
            do.call("cbind",
                    parallel::mclapply(split_index(B, ncpus),
                                       FUN = montecarlo, mc.cores = ncpus,
                                       mc.allow.recursive = FALSE))
        } else {
            if (is.null(cl)) {
                ## has a default cluster been registered?
                ## see parallel:::defaultCluster
                cl <- get("default",
                          envir = get(".reg", envir = getNamespace("parallel"),
                                      inherits = FALSE),
                          inherits = FALSE)
                if (is.null(cl)) {
                    ## no default cluster, so setup a PSOCK cluster
                    cl <- parallel::makePSOCKcluster(ncpus)
                    on.exit(parallel::stopCluster(cl), add = TRUE) # clean-up
                }
            }
            if (RNGkind()[1L] == "L'Ecuyer-CMRG")
                ## distribute streams (using master process) for reproducibility
                parallel::clusterSetRNGStream(cl)
            ncpus <- as.integer(length(cl))
            if (ncpus < 2L)
                warning("parallel operation requires at least two processes")
###            Bp <- split_index(B, ncpus) # distribute workload evenly
###            RET <- parallel::clusterApply(cl, x = Bp, fun = montecarlo)
            do.call("cbind",
                    parallel::clusterApply(cl, x = split_index(B, ncpus),
                                           fun = montecarlo))
        }
###        do.call("cbind", RET)
    }
}
