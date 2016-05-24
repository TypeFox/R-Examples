summary.plsRmodel <- function(object, ...)
{
res <- list(call=object$call)
class(res) <- "summary.plsRmodel"
res
}
