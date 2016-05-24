predict.hdglm <-
function(object, newdata, se.fit = FALSE, scale = NULL, df = Inf,
     interval = c("none", "confidence", "prediction"),
     level = .95,  type = c("response", "terms"),
     terms = NULL, na.action = na.pass, weights = 1, ...)
{
    tt <- object
    # Simple / less complete than lm version:
    if(!inherits(object, "hdlm"))
        warning("calling predict.hdlm(<fake-lm-object>) ...")
    if(missing(newdata) || is.null(newdata)) {
        X <- model.matrix(object)
    } else {
        X <- model.matrix(object$fitted.values ~ newdata)
        if(variable.names(object)[[1]] != "(Intercept)") X <- X[,-1]
    }

    beta <- coef(object)
    ynew <- X %*% beta

    if(object$family == 'binomial') {
        ynew <- 1 / (1 + exp(ynew))
    } else if (object$family == 'poisson') {
        ynew <- exp(ynew)
    }

    return(ynew)
}

