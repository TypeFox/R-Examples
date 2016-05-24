summary.plsRbetamodel <- function(object, ...)
{
res <- list(call=object$call)
class(res) <- "summary.plsRbetamodel"
res
}