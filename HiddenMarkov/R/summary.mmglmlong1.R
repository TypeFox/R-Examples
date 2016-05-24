"summary.mmglmlong1" <-
function (object, ...) 
{
    N <- length(table(object$longitude))
    n <- length(object$y)/N
    list(delta=object$delta, Pi=object$Pi, nonstat=object$nonstat,
         beta=object$beta,
         sigma=object$sigma, glmfamily=object$glmfamily,
         n=n, N=N)
}

