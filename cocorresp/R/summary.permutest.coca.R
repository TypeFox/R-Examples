"summary.permutest.coca" <- function(object, ...)
{
    retval <- list(pval = object$pval, permstat = object$permstat,
                   inertia = object$inertia,
                   fitax = object$fitax, pcent.fit = object$pcent.fit,
                   n.axes = object$n.axes,
                   total.inertia = object$total.inertia,
                   call = object$call)
    ## Ychi1 = object$Ychi1, Ychi2 = object$Ychi2)
    class(retval) <- "summary.permutest.coca"
    retval
}

