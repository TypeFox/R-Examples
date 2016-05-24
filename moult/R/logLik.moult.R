
logLik.moult <- function (object, ...)
{ structure(object$loglik, df = object$n - object$df.residual,
            class = "logLik")
}
