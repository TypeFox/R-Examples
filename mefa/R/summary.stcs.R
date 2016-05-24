`summary.stcs` <-
function(object, ...)
{
out <- summary.mefa(mefa(object, ...))
out$call <- attr(out, "call")
class(out) <- c("summary.stcs", "summary.mefa")
return(out)
}

