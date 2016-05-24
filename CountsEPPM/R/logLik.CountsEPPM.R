logLik.CountsEPPM <-
function (object, ...) {
      structure(object$loglikelihood, df = nrow(object$estses),
                nobs = nrow(object$covariates.matrix.mean),
                class = "logLik") }
