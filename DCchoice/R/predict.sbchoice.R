predict.sbchoice <- function(object, newdata = NULL,
    type = c("utility", "probability"), ...)
{

    type <- match.arg(type)

    COEF <- matrix(object$coefficients, ncol = 1)

    if(is.null(newdata)) {
        X <- model.matrix(object$formula, data = object$data.name, rhs = 1:2)
        V <- X %*% COEF
    } else {
        formula <- delete.response(object$terms)
        mf.newX <- model.frame(formula, newdata, xlev = object$xlevels)
        mm.newX <- model.matrix(formula, mf.newX, contrasts.arg = object$contrasts)
        V <- mm.newX %*% COEF
    }

    dist <- object$distribution
    if(dist == "logistic" | dist == "log-logistic") {
        P <- plogis(-V, lower.tail = FALSE, log.p = FALSE)
    } else if (dist == "normal" | dist == "log-normal") {
        P <- pnorm(-V, lower.tail = FALSE, log.p = FALSE)
    } else if (dist == "weibull") {
        P <- pweibull(exp(-V), shape = 1, lower.tail = FALSE, log.p = FALSE)
    }

    if(type == "utility") {
        return(as.vector(V))
    } else {
        return(as.vector(P))
    }
}
