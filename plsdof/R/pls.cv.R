pls.cv<-function (X, y, k = 10, groups = NULL, m = ncol(X), use.kernel = FALSE, 
    compute.covariance = FALSE,method.cor="pearson") 
{
    n <- nrow(X)
    p <- ncol(X)
    if (is.null(groups) == FALSE) {
        f = as.factor(groups)
        k = length(levels(f))
        my.names = levels(f)
        all.folds <- split(1:n, f)
    }
    if (is.null(groups) == TRUE) {
        f <- rep(1:k, length = n)
        my.names <- 1:k
       all.folds <- split(sample(1:n), f)
    }
    
    ntrain = vector(length = k)
    for (i in 1:k) {
        ntrain[i] = n - length(all.folds[[i]])
    }
    ntrain.min = min(ntrain)
    m = min(m, ntrain.min - 1, p)
    cv.error.matrix = matrix(0, k, m + 1)
    rownames(cv.error.matrix) = my.names
    colnames(cv.error.matrix) = 0:m
    cor.error.matrix<-cv.error.matrix
    for (i in seq(k)) {
        omit <- all.folds[[i]]
        Xtrain = X[-omit, , drop = FALSE]
        ytrain = y[-omit]
        Xtest = X[omit, , drop = FALSE]
        ytest = y[omit]
        pls.object <- pls.model(Xtrain, ytrain, m = m, Xtest = Xtest, 
            ytest = ytest, compute.DoF = FALSE, use.kernel = use.kernel,method.cor=method.cor)
        cv.error.matrix[i, ] <- pls.object$mse
        cor.error.matrix[i, ] <- pls.object$cor
    }
    cv.error = apply(cv.error.matrix, 2, mean)
    cor.error<-apply(cor.error.matrix,2,mean)
    m.opt <- which.min(cv.error) - 1
    m.opt.cor<-which.max(cor.error) - 1
    if (compute.covariance == TRUE) {
        use.kernel = FALSE
    }
    pls.object <- pls.model(X, y, m = max(m.opt, m.opt.cor,1), use.kernel = use.kernel, 
        compute.DoF = compute.covariance, compute.jacobian = compute.covariance)
    intercept <- pls.object$intercept[m.opt + 1]
    coefficients <- pls.object$coefficients[, m.opt + 1]
    covariance <- pls.object$covariance
    intercept.cor <- pls.object$intercept[m.opt.cor + 1]
    coefficients.cor <- pls.object$coefficients[, m.opt.cor + 1]
    if (compute.covariance == TRUE) {
        #covariancve.cor<-covariance[m.opt.cor + 1, , ]
        covariance <- covariance[m.opt + 1, , ]
    }
    outlist = list(cv.error.matrix = cv.error.matrix, cor.error.matrix=cor.error.matrix,cv.error = cv.error, cor.error=cor.error,
        m.opt = m.opt, m.opt.cor=m.opt.cor,covariance = covariance, intercept = intercept, intercept.cor=intercept.cor,
        coefficients = coefficients,coefficients.cor=coefficients.cor)
    class(outlist) = "plsdof"
    return(outlist)
}


