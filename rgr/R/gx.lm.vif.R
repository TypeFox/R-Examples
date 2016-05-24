gx.lm.vif <-
function (object, ...) 
{
    V <- summary(object)$cov.unscaled
    Vi <- crossprod(model.matrix(object))
    nam <- names(coef(object))
    if (k <- match("(Intercept)", nam, nomatch = FALSE)) {
        v1 <- diag(V)[-k]
        v2 <- (diag(Vi)[-k] - Vi[k, -k]^2/Vi[k, k])
        nam <- nam[-k]
    }
    else {
        v1 <- diag(V)
        v2 <- diag(Vi)
        warning(" No intercept term detected - Results may surprise.")
    }
    structure(v1 * v2, names = nam)
}
