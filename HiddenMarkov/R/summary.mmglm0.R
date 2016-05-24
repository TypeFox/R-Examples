"summary.mmglm0" <-
function (object, ...) 
{
    list(variable.names=names(object$x),
         delta=object$delta, Pi=object$Pi, nonstat=object$nonstat,
         beta=object$beta,
         sigma=object$sigma, family=object$family,
         glmformula=object$glmformula,
         link=object$link, n=length(object$x$y))
}

