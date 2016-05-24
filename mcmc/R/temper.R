
temper <- function(obj, initial, neighbors, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, debug = FALSE, parallel = FALSE, ...)
UseMethod("temper")

temper.tempering <- function(obj, initial, neighbors, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, debug = FALSE, parallel = FALSE, ...)
{
    if (missing(initial)) initial <- obj$final
    if (missing(neighbors)) neighbors <- obj$neighbors
    if (missing(nbatch)) nbatch <- obj$nbatch
    if (missing(blen)) blen <- obj$blen
    if (missing(nspac)) nspac <- obj$nspac
    if (missing(scale)) scale <- obj$scale
    if (missing(debug)) debug <- obj$debug
    if (missing(parallel)) parallel <- obj$parallel
    assign(".Random.seed", obj$final.seed, .GlobalEnv)
    if (missing(outfun)) {
        if (is.null(obj$outfun)) {
            temper.function(obj$lud, initial, neighbors, nbatch, blen,
                nspac, scale, debug = debug, parallel = parallel, ...)
        } else {
            temper.function(obj$lud, initial, neighbors, nbatch, blen,
                nspac, scale, obj$outfun, debug = debug, parallel = parallel,
                ...)
        }
    } else {
        temper.function(obj$lud, initial, neighbors, nbatch, blen,
            nspac, scale, outfun, debug = debug, parallel = parallel, ...)
    }
}

temper.function <- function(obj, initial, neighbors, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, debug = FALSE, parallel = FALSE, ...)
{
    if (! exists(".Random.seed")) runif(1)
    saveseed <- .Random.seed
    func1 <- function(state) obj(state, ...)
    env1 <- environment(fun = func1)
    if (missing(outfun)) {
        func2 <- NULL
        env2 <- NULL
        outfun <- NULL
    } else if (is.function(outfun)) {
        func2 <- function(state) outfun(state, ...)
        env2 <- environment(fun = func2)
    }

    stopifnot(is.numeric(initial))
    storage.mode(initial) <- "double"

    if (is.list(scale)) {
        for (i in 1:length(scale)) {
            stopifnot(is.numeric(scale[[i]]))
            storage.mode(scale[[i]]) <- "double"
        }
    } else {
        stopifnot(is.numeric(scale))
        storage.mode(scale) <- "double"
    }

    out.time <- system.time(
    out <- .Call("temper", func1, initial, neighbors, nbatch, blen, nspac,
        scale, func2, debug, parallel, env1, env2, PACKAGE = "mcmc")
    )
    result <- structure(c(list(lud = obj, initial = initial,
        neighbors = neighbors, nbatch = nbatch, blen = blen, nspac = nspac,
        scale = scale, outfun = outfun, debug = debug, parallel = parallel,
        initial.seed = saveseed, final.seed = .Random.seed, time = out.time),
        out), class = c("mcmc", "tempering"))
    return(result)
}

