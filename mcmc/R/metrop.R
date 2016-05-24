
metrop <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, debug = FALSE, ...)
UseMethod("metrop")

metrop.metropolis <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, debug = FALSE, ...)
{
    if (missing(nbatch)) nbatch <- obj$nbatch
    if (missing(blen)) blen <- obj$blen
    if (missing(nspac)) nspac <- obj$nspac
    if (missing(scale)) scale <- obj$scale
    if (missing(debug)) debug <- obj$debug
    assign(".Random.seed", obj$final.seed, .GlobalEnv)
    if (missing(outfun)) {
        if (is.null(obj$outfun)) {
            metrop.function(obj$lud, obj$final, nbatch, blen,
                nspac, scale, debug = debug, ...)
        } else {
            metrop.function(obj$lud, obj$final, nbatch, blen,
                nspac, scale, obj$outfun, debug, ...)
        }
    } else {
        metrop.function(obj$lud, obj$final, nbatch, blen,
            nspac, scale, outfun, debug, ...)
    }
}

metrop.function <- function(obj, initial, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, debug = FALSE, ...)
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
    } else {
        func2 <- outfun
        env2 <- NULL
    }

    out.time <- system.time(
    out <- .Call("metrop", func1, initial, nbatch, blen, nspac,
        scale, func2, debug, env1, env2, PACKAGE = "mcmc")
    )
    out$initial.seed <- saveseed
    out$final.seed <- .Random.seed
    out$time <- out.time
    out$lud <- obj
    out$nbatch <- nbatch
    out$blen <- blen
    out$nspac <- nspac
    out$scale <- scale
    out$outfun <- outfun
    out$batch <- t(out$batch)
    out$debug <- debug
    if (! is.null(out$current)) out$current <- t(out$current)
    if (! is.null(out$proposal)) out$proposal <- t(out$proposal)
    if (! is.null(out$z)) out$z <- t(out$z)
    class(out) <- c("mcmc", "metropolis")
    return(out)
}

