"as.dthmm" <-
function (object)
{
    #   turns mmglm0 object into dthmm
    if (object$family=="binomial" | object$family=="poisson") discrete <- TRUE
    else discrete <- FALSE
    x <- list(x=object$x$y, Pi=object$Pi, delta=object$delta,
              distn=c("glm", object$family, object$link),
              pm=list(beta0=object$beta[1,], beta1=object$beta[2,], sigma=object$sigma),
              pn=list(x1=object$x$x1, size=object$x$size), nonstat=object$nonstat,
              discrete=discrete)
    x$glmdata <- object$x
    x$glmformula <- object$glmformula
    class(x) <- c("dthmm")
    return(x)
}

