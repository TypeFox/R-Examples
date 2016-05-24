"summary.mmglm1" <-
function (object, ...) 
{
    list(delta=object$delta, Pi=object$Pi, nonstat=object$nonstat,
         beta=object$beta,
         sigma=object$sigma, glmfamily=object$glmfamily,
         n=length(object$y))
}

