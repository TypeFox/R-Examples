"summary.FitARMA" <-
function (object, ...) 
{
    LL <- object$loglikelihood
    k <- object$order[1]+object$order[3]
    if (!is.null(object$demean) && object$demean) 
        k <- k + 1
    n <- length(object$res)
    aic <- -2 * LL + 2 * k
    bic <- -2 * LL + log(n) * k
    dati <- object$DataTitle
    if (!is.null(dati)) 
        cat(dati, fill = TRUE)
    modti <- object$ModelTitle
    if (object$MeanMLE) 
        modti <- paste(modti, " With mean MLE.")
    cat(modti, fill = TRUE)
    cat(paste("length of series =", n, ",  number of parameters =", 
        k), fill = TRUE)
    cat(paste("loglikelihood =", round(LL, 2), ",  aic =", round(aic, 
        1), ",  bic = ", round(bic, 1)), fill = TRUE)
}

