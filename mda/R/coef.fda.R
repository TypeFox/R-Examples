coef.fda <-
function (object, ...) 
{
    fit <- object$fit
    Coefs <- fit$coef
    if (is.null(Coefs)) 
        stop("No explicit coefficients in this formulation")
    dimension <- object$dimension
    Coefs <- Coefs %*% object$theta[, seq(dimension), drop = FALSE]
    lambda <- object$values
    alpha <- sqrt(lambda[seq(dimension)])
    sqima <- sqrt(1 - lambda[seq(dimension)])
    Coefs <- scale(Coefs, FALSE, sqima * alpha)
    attr(Coefs,"scaled:scale")=NULL
    Coefs
}

