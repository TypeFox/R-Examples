logLik.SemiParSampleSel <- function (object, ...){

    if (length(list(...))) warning("extra arguments discarded")
    
    lk  <- object$logLik
    
    attr(lk, "nobs") <- object$n
    attr(lk, "df") <- object$t.edf
    class(lk) <- "logLik"
    lk
    
}
