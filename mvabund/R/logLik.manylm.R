################################################################################
# Extract the Log-Likelihood of a manylm object                                #
################################################################################

logLik.manylm <- 
function (object, REML = FALSE, ...) {
    res <- object$residuals
    p <- object$rank
    N <- NROW(res)
    q <- NCOL(res)
    if (is.null(w <- object$weights)) {
        w <- rep.int(1, N)
    } else {
        excl <- w == 0
        if (any(excl)) {
            res <- res[!excl,]
            N <- NROW(res)
            w <- w[!excl]
        }
    }
    N0 <- N
    if (REML) 
        N <- N - p
        
    val <- 0.5 * (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) + 
        log(colSums(w * res^2))) )
    if (REML) 
        val <- val - sum(log(abs(diag(object$qr$qr)[1:p])))
        
    attr(val, "nall") <- N0
    attr(val, "nobs") <- N
    attr(val, "vars") <- q
    attr(val, "df")   <- p + 1
    class(val) <- "logLik"
    val

}
