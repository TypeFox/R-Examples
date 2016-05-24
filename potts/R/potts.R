
potts <- function(obj, param, nbatch, blen = 1, nspac = 1,
    boundary = c("torus", "free", "condition"), debug = FALSE,
    outfun = NULL, ...)
UseMethod("potts")

potts.potts <- function(obj, param, nbatch, blen = 1, nspac = 1,
    boundary = c("torus", "free", "condition"), debug = FALSE,
    outfun = NULL, ...)
{
    boundary <- match.arg(boundary)
    if (missing(param)) param <- obj$param
    if (missing(nbatch)) nbatch <- obj$nbatch
    if (missing(blen)) blen <- obj$blen
    if (missing(nspac)) nspac <- obj$nspac
    if (missing(boundary)) boundary <- obj$boundary
    if (missing(debug)) debug <- obj$debug
    if (missing(outfun)) outfun <- obj$outfun
    initial <- obj$final
    assign(".Random.seed", obj$final.seed, .GlobalEnv)
    potts.raw(initial, param, nbatch, blen, nspac, boundary, debug, outfun, ...)
}

potts.raw <- function(obj, param, nbatch, blen = 1, nspac = 1,
    boundary = c("torus", "free", "condition"), debug = FALSE,
    outfun = NULL, ...)
{
    boundary <- match.arg(boundary)

    if (! exists(".Random.seed")) runif(1)
    saveseed <- .Random.seed

    initial.info <- inspectPotts(obj)
    stopifnot(is.numeric(param))
    stopifnot(all(is.finite(param)))
    if(length(param) != initial.info$ncolor + 1)
        stop("length(param) not number of colors + 1")
    nrow <- initial.info$nrow
    ncol <- initial.info$ncol

    stopifnot(is.numeric(nbatch))
    stopifnot(length(nbatch) == 1)
    stopifnot(nbatch == as.integer(nbatch))
    nbatch <- as.integer(nbatch)
    stopifnot(nbatch > 0)

    stopifnot(is.numeric(blen))
    stopifnot(length(blen) == 1)
    stopifnot(blen == as.integer(blen))
    blen <- as.integer(blen)
    stopifnot(blen > 0)

    stopifnot(is.numeric(nspac))
    stopifnot(length(nspac) == 1)
    stopifnot(nspac == as.integer(nspac))
    nspac <- as.integer(nspac)
    stopifnot(nspac > 0)

    niter <- nbatch * blen * nspac

    boundary.code <- match(boundary, c("torus", "free", "condition"))

    stopifnot(is.logical(debug))
    stopifnot(length(debug) == 1)
    if (debug) {
        pstate <- array(as.integer(0), dim = c(niter, nrow, ncol))
        hstate <- array(as.integer(0), dim = c(niter, nrow, ncol))
        vstate <- array(as.integer(0), dim = c(niter, nrow, ncol))
        patch <- array(as.integer(0), dim = c(niter, nrow, ncol))
        hunif <- array(as.double(-1), dim = c(niter, nrow, ncol))
        vunif <- array(as.double(-1), dim = c(niter, nrow, ncol))
        punif <- matrix(as.double(-1), niter, nrow * ncol)
    } else {
        pstate <- array(as.integer(0), dim = c(1, 1, 1))
        hstate <- array(as.integer(0), dim = c(1, 1, 1))
        vstate <- array(as.integer(0), dim = c(1, 1, 1))
        patch <- array(as.integer(0), dim = c(1, 1, 1))
        hunif <- array(as.double(-1), dim = c(1, 1, 1))
        vunif <- array(as.double(-1), dim = c(1, 1, 1))
        punif <- array(as.double(-1), dim = c(1, 1, 1))
        punif <- matrix(as.double(-1), 1, 1)
    }

    if (is.null(outfun)) {
        .C("outfun_shutdown", PACKAGE = "potts")
        nout <- length(param)
    } else {
        func2 <- function(tt) outfun(tt, ...)
        env2 <- environment(fun = func2)
        .Call("outfun_setup", func2, env2)
        nout <- .C("outfun_len_init", x = obj,
            code = as.integer(boundary.code), nout = integer(1),
            PACKAGE = "potts")$nout
    }

    out.time <- system.time(
    out <- .C("potts", final = obj, param = as.double(param),
        nbatch = nbatch, blen = blen, nspac = nspac,
        code = as.integer(boundary.code),
        batch = matrix(as.double(0), nrow = as.integer(nout),
            ncol = as.integer(nbatch)),
        debug = debug, pstate = pstate, hstate = hstate, vstate = vstate,
        patch = patch, hunif = hunif, vunif = vunif, punif = punif,
        PACKAGE = "potts")
    )
    .C("outfun_shutdown", PACKAGE = "potts")

    if (debug) {
        return(structure(list(initial.seed = saveseed,
        final.seed = .Random.seed, initial = obj, final = out$final,
        param = param, nbatch = nbatch, blen = blen, nspac = nspac,
        boundary = boundary, batch = t(out$batch), time = out.time,
        debug = TRUE, pstate = out$pstate, hstate = out$hstate,
        vstate = out$vstate, patch = out$patch, hunif = out$hunif,
        vunif = out$vunif, punif = out$punif), class = "potts"))
    }

    return(structure(list(initial.seed = saveseed, final.seed = .Random.seed,
        initial = obj, final = out$final, param = param, nbatch = nbatch,
        blen = blen, nspac = nspac, boundary = boundary,
        batch = t(out$batch), time = out.time, debug = FALSE),
        class = "potts"))
}

