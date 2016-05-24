pool.RR <-
function (object, method = "plain")
{
    call <- match.call()
#    if (!is.mira(object))
#        stop("The object must have class 'mira'")
    if ((m <- length(object)) < 2)
        stop("At least two imputations are needed for pooling.\n")
    analyses <- object
    mess <- try(vcov(analyses[[1]]), silent = TRUE)
    if (inherits(mess, "try-error"))
        stop("Object has no vcov() function.")
    k <- length(coef(analyses[[1]]))
    names <- names(coef(analyses[[1]]))

    qhat <- matrix(NA, nrow = m, ncol = k, dimnames = list(1:m,
        names))
    u <- array(NA, dim = c(m, k, k), dimnames = list(1:m, names,
        names))
    for (i in 1:m) {
        fit <- analyses[[i]]
            if (class(fit)[1] == "lme") {
                qhat[i, ] <- fit$coefficients$fixed
            }
            else qhat[i, ] <- coef(fit)
            u[i, , ] <- vcov(fit)
    }
    qbar <- apply(qhat, 2, mean)
    ubar <- apply(u, c(2, 3), mean)
    e <- qhat - matrix(qbar, nrow = m, ncol = k, byrow = TRUE)
    b <- (t(e) %*% e)/(m - 1)
    t <- ubar + (1 + 1/m) * b
    r <- (1 + 1/m) * diag(b/ubar)
    f <- (1 + 1/m) * diag(b/t)
    df <- (m - 1) * (1 + 1/r)^2
    lambda <- (1 + 1/m) * diag(b/t)
    if (method == "smallsample") {
        cls <- class(fit)[1]
        if (cls == "lme")
            dfc <- fit$fixDF[["X"]]
        else if (cls == "multinom")
            dfc <- fit$edf
        else if (is.null(fit$df.residual))
            stop("Cannot extract df from object of class ", cls,
                ". Use pool(.., method=\"plain\").")
        else dfc <- fit$df.residual
        df <- dfc/((1 - (f/(m + 1)))/(1 - f) + dfc/df)
    }
    names(r) <- names(df) <- names(f) <- names
    fit <- list(call = call, call1 = object[[1]]$call, call2 = NA,
        nmis = NA, m = m, qhat = qhat, u = u, qbar = qbar,
        ubar = ubar, b = b, t = t, r = r, df = df, fmi = f, lambda=lambda)
    oldClass(fit) <- c("mipo", oldClass(object))
    return(fit)
}

