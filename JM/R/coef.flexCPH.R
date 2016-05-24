coef.flexCPH <-
function (object, include.splineCoefs = FALSE, ...) {
    if (include.splineCoefs) unlist(object$coefficients) else object$coefficients$betas
}
