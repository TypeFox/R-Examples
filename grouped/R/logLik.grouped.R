"logLik.grouped" <-
function(object, ...){
    dets <- object$details 
    out <- dets$logLik
    attributes(out) <- list(nall = dets$n, nobs = dets$n, "df" = length(object$coef), class = "logLik")
    out
}

