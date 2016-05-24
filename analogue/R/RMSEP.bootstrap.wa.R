RMSEP.bootstrap.wa <- function(object,
                               type = c("birks1990", "standard"),
                               ...) {
    if (missing(type))
        type <- "birks1990"
    type <- match.arg(type)
    if (type == "birks1990")
        rmsep <- object$performance$rmsep
    else rmsep <- sqrt(mean(object$model$resid^2))
    return(rmsep)
}
