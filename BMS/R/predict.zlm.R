predict.zlm <-
function (object, newdata = NULL, se.fit = FALSE, ...) 
{
    if (!is(object, "zlm")) {
        stop("you need to provide a zlm object")
        return()
    }
    betas = object$coefficients[-1, drop = FALSE]
    alpha = object$coefficients[[1]]
    if (is.null(newdata)) {
        newX <- as.matrix(object$model[, -1, drop = FALSE])
    }
    else {
        newX = as.matrix(newdata)
        if (!is.numeric(newX)) 
            stop("newdata must be numeric!")
        if (is.vector(newdata)) 
            newX = matrix(newdata, 1)
        if (ncol(newX) != length(betas)) {
            if (ncol(newX) == length(betas) + 1) {
                newX = newX[, -1, drop = FALSE]
            }
            else {
                stop("newdata must be a matrix or data.frame with ", 
                  length(betas), " columns.")
            }
        }
    }
    if (!se.fit) 
        return(as.vector(newX %*% betas) + alpha)
    yXdata <- as.matrix(object$model)
    oldXraw <- yXdata[, -1, drop = FALSE]
    if (!is.null(colnames(newX)) && !is.null(colnames(oldXraw))) {
        if (all(colnames(oldXraw) %in% colnames(newX)) && !all(colnames(oldXraw) == 
            colnames(newX))) {
            warning("argument newdata had to be reordered according to its column names. Consider submitting the columns of newdata in the right order.")
            newX = newX[, colnames(oldXraw), drop = FALSE]
        }
    }
    yraw <- yXdata[, 1, drop = TRUE]
    N <- length(yraw)
    k <- ncol(oldXraw)
    oldXmeans <- colMeans(oldXraw)
    oldXdm <- oldXraw - matrix(oldXmeans, N, k, byrow = TRUE)
    newXdm <- newX - matrix(oldXmeans, nrow(newX), k, byrow = TRUE)
    xtxinv = chol2inv(chol(crossprod(oldXdm)))
    xtxinv_xf = tcrossprod(xtxinv, newXdm)
    xf_xx_xf = unlist(lapply(1:nrow(newXdm), function(x) {
        crossprod(newXdm[x, ], xtxinv_xf[, x])[[1L]]
    }))
    bvar = object$coef2moments[-1] - object$coefficients[-1]^2
    bvar_factor = bvar[[1L]]/xtxinv[[1L]]
    yty = as.vector(crossprod(yraw) - N * mean(yraw)^2)
    r2 = 1 - object$olsres$ymy/yty
    if (object$gprior.info$gtype == "hyper") {
        f21a = object$gprior.info$hyper.parameter
        f21_recover = exp((object$marg.lik) + (N - 1)/2 * log(yty) + 
            log((k + f21a - 2)/(f21a - 2)))
        res_scale = yty/(N - 3)/(N - 1 - k - f21a) * ((N - 3) * 
            (1 - r2) - (k + f21a - 2)/f21_recover)
        svar_woscale = res_scale/N + bvar_factor * xf_xx_xf
    }
    else {
        sf = object$gprior.info$shrinkage.moments[[1]]
        res_scale = (1 - sf * r2) * yty/(N - 3)
        svar_woscale = res_scale/N + bvar_factor * xf_xx_xf
    }
    reslist = list()
    reslist$fit <- as.vector(newX %*% betas) + alpha
    reslist$std.err <- sqrt(svar_woscale + res_scale)
    reslist$se.fit <- sqrt(svar_woscale)
    reslist$residual.scale <- sqrt(res_scale)
    return(reslist)
}
