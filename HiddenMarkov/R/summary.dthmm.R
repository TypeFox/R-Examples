"summary.dthmm" <-
function (object, ...) 
{
    list(delta=object$delta, Pi=object$Pi, nonstat=object$nonstat,
         distn=object$distn, pm=object$pm, discrete=object$discrete,
         n=length(object$x))
}

