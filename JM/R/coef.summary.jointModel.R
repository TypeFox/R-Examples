coef.summary.jointModel <-
function (object, ...) {
    if (!inherits(object, "summary.jointModel"))
        stop("Use only with 'summary.jointModel' objects.\n")
    coefsY <- object$'CoefTable-Long'
    coefsT <- object$'CoefTable-Event'
    list("Longitudinal" = coefsY, "Event" = coefsT)
}
