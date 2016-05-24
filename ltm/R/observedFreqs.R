observedFreqs <-
function (object, Y) {
    if (!class(object) %in% c("grm", "gpcm", "ltm", "rasch", "tpm"))
        stop("'object' must inherit from either class 'grm', class 'gpcm', class 'ltm', class 'rasch' or class 'tpm'.\n")
    X <- object$patterns$X
    Obs <- object$patterns$obs
    patsX <- apply(X, 1, paste, collapse = "/")
    patsY <- apply(Y, 1, paste, collapse = "/")
    obs <- numeric(nrow(Y))
    ind <- match(patsY, patsX)
    na <- !is.na(ind)
    obs[na] <- Obs[ind[na]]
    obs
}
