raw.pvalues <-
function (data, x = 1, test = x[1]) {
    if (class(data)[1] != "FAMTdata")
        stop("Class of data should be FAMTdata")
    nbcovar = ncol(data$covariates)
    if (any(!is.element(x, 1:nbcovar)))
        stop(paste("x should be a subset of 1:", nbcovar, sep = ""))
    if (any(!is.element(test, 1:nbcovar)))
        stop(paste("test should be one of ", x, sep = ""))
    n = nrow(data$covariates)
    m = nrow(data$expression)
    y = t(data$expression)
    xtest = which(is.element(x,test))
    if (all(data$covariates[, 2] == rep(1, n))) {
        m = apply(data$expression, 1, mean)
        s = apply(data$expression, 1, sd)
        Ftest = (m/s)^2 * n
        dfr1 = n - 1
        dfr0 = n
        pval = pf(Ftest, df1 = dfr0 - dfr1, df2 = dfr1, lower.tail = FALSE)
    }
    if (any(data$covariates[, 2] != rep(1, n))) {
        covar = data$covariates[, c(x[-xtest], test)]
        X1 = model.matrix(y ~ ., data = data.frame(y = rep(1,
            n), covar))
        P1 = diag(n) - X1 %*% solve(t(X1) %*% X1) %*% t(X1)
        SCER1 = apply(y * (P1 %*% y), 2, sum)
        dfr1 = sum(diag(P1))
        X0 = matrix(rep(1, n), ncol = 1)
        if (length(x) > 1)
            X0 = model.matrix(y ~ ., data = data.frame(y = rep(1,
                n), covar[, -ncol(covar)]))
        P0 = diag(n) - X0 %*% solve(t(X0) %*% X0) %*% t(X0)
        SCER0 = apply(y * (P0 %*% y), 2, sum)
        dfr0 = sum(diag(P0))
        Ftest = ((SCER0 - SCER1)/(dfr0 - dfr1))/(SCER1/dfr1)
        pval = pf(Ftest, df1 = dfr0 - dfr1, df2 = dfr1, lower.tail = FALSE)
    }
    list(pval = pval, test = Ftest, resdf = dfr1)
}
