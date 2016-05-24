residuals.BTm <- function(object, type = c("deviance", "pearson", "working",
                          "response", "partial", "grouped"), by = object$id,
                          ...) {
    type <- match.arg(type)
    if (type != "grouped") return(NextMethod())

    ## for glm, lm would just be
    ## X <- model.matrix(formula, data = object$data)
    formula <- as.formula(paste("~", by, "- 1"))
    mt <- terms(formula)
    mf1 <- model.frame(mt, data = c(object$player1, object$data))
    X1 <- model.matrix(mt, data = mf1)
    mf2 <- model.frame(mt, data = c(object$player2, object$data))
    X2 <- model.matrix(mt, data = mf2)
    X <- X1 - X2

    r <- object$residuals  ## the "working" residuals
    w <- object$weights
    total.resid <- crossprod(X, r * w)
    total.weight <- crossprod(abs(X), w)
    result <- total.resid / total.weight
    attr(result, "weights") <- total.weight
    result
}
