simulate.gpcm <-
function (object, nsim = 1, seed = NULL, ...) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    if (is.null(seed)) 
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    n <- nrow(object$X)
    ans <- lapply(seq_len(nsim), function (i) {
        rmvordlogis(nrow(object$X), object$coefficients, object$IRT.param, model = "gpcm")
    })
    attr(ans, "seed") <- RNGstate
    ans
}
