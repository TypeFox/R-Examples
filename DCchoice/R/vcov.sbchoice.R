vcov.sbchoice <- function(object, ...)
{
    if (object$dist == "weibull") {
        solve(object$glm.out$hessian)
    } else {
        summary.glm(object$glm.out)$cov.scaled
    }
}
