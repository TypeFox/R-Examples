logLik.modTempEff<-function (object, ...){
    if (length(list(...))) warning("extra arguments discarded")
    p <- sum(object$edf)
    val <- p - object$aic/2
    attr(val, "df") <- p
    class(val) <- "logLik"
    val
}