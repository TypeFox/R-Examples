"summary.mmpp" <-
function (object, ...) 
{
    list(delta=object$delta, Q=object$Q, nonstat=object$nonstat,
         lambda=object$lambda, n=length(object$tau))
}

