summary.plsRcoxmodel <- function(object, ...)
{
res <- list(call=object$call)
class(res) <- "summary.plsRcoxmodel"
res
}