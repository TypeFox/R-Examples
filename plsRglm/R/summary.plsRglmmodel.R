summary.plsRglmmodel <- function(object, ...)
{
res <- list(call=object$call)
class(res) <- "summary.plsRglmmodel"
res
}