residualsFAMT <-
function (data, x = 1, test = x, pvalues = NULL) {
    if (class(data)[1] != "FAMTdata")
        stop("Class of data should be FAMTdata")
    nbcovar = ncol(data$covariates)
    if (any(!is.element(x, 1:nbcovar)))
        stop(paste("x should be a subset of 1:", nbcovar, sep = ""))
    if (any(!is.element(test, 1:nbcovar)))
        stop(paste("test should be one of ", x, sep = ""))
    threshold = 0.05
    n = ncol(data$expression)
    m = nrow(data$expression)
    if (length(pvalues) == 0)
        pvalues = raw.pvalues(data, x, test)$pval
    SelectH0 = (1:m)[pvalues > threshold]
    y = t(data$expression)
    xtest = which(is.element(x,test))
    covar = data$covariates[, c(x[-xtest], test)]
    if (all(data$covariates[, 2] == rep(1, n))) {
        P0 = diag(n)
        P1 = diag(n) - (1/n) * matrix(1, n, n)
    }
    if (any(data$covariates[, 2] != rep(1, n))) {
        X0 = matrix(rep(1, n), ncol = 1)
        if (length(x) > 1)
            X0 = model.matrix(y ~ ., data = data.frame(y = rep(1,
                n), covar[, -ncol(covar)]))
        P0 = diag(n) - X0 %*% solve(t(X0) %*% X0) %*% t(X0)
        X1 = model.matrix(y ~ ., data = data.frame(y = rep(1,
            n), covar))
        P1 = diag(n) - X1 %*% solve(t(X1) %*% X1) %*% t(X1)
    }
    res = P0 %*% y
    res[, -SelectH0] = P1 %*% y[, -SelectH0]
    list(residuals = res, SelectH0 = SelectH0, P = P1)
}
