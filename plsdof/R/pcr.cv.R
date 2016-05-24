 pcr.cv<-function (X, y, k = 10, m = min(ncol(X), nrow(X) - 1), groups = NULL, 
    scale = TRUE, eps = 1e-06, plot.it = FALSE, compute.jackknife = TRUE,method.cor="pearson",supervised=FALSE) 
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
    m = min(m, ntrain.min - 1)
    cv.error.matrix = matrix(0, k, m + 1)
    rownames(cv.error.matrix) = my.names
    colnames(cv.error.matrix) = 0:m
    cor.error.matrix<-cv.error.matrix
    coefficients.jackknife = NULL
    if (compute.jackknife == TRUE) {
        coefficients.jackknife <- array(dim = c(p, m + 1, k))
    }
    for (i in seq(k)) {
        omit <- all.folds[[i]]
        Xtrain = X[-omit, , drop = FALSE]
        ytrain = y[-omit]
        Xtest = X[omit, , drop = FALSE]
        ytest = y[omit]
        pcr.object <- pcr(Xtrain, ytrain, scale, m, eps,supervised)
        res <- matrix(, length(ytest), m + 1)
        if (compute.jackknife == TRUE) {
            coefficients.jackknife[, , i] = pcr.object$coefficients
        }
        for (j in 1:(m + 1)) {
            ytest.j<-pcr.object$intercept[j]+ Xtest%*%pcr.object$coefficients[, j]
            res[, j] <- ytest - ytest.j
            if (j>1){
            cor.error.matrix[i,j]<-cor(ytest,ytest.j,method=method.cor)
            }
        }
        cv.error.matrix[i, ] <- apply(res^2, 2, mean)
    }
 #   cat(paste("cv complete \n"))
    cv.error <- apply(cv.error.matrix, 2, mean)
    cor.error <- apply(cor.error.matrix, 2, mean)
    m.opt <- which.min(cv.error) - 1
 #  cat(paste("mse ",m.opt,"\n"))
    m.opt.cor <- which.max(cor.error) - 1
#	cat(paste("mse ",m.opt,"\n"))
    if (plot.it == TRUE) {
        plot(0:m, cv.error, type = "l")
    }
    m.min<-max(2,m.opt,m.opt.cor)
    pcr.object <- pcr(X, y, scale=scale,m=m.min, eps = eps,supervised=supervised)
    coefficients <- pcr.object$coefficients[, m.opt + 1]
    intercept.cor <- pcr.object$intercept[m.opt.cor + 1]
    coefficients.cor <- pcr.object$coefficients[, m.opt.cor + 1]
    intercept <- pcr.object$intercept[m.opt + 1]
    return(list(intercept = intercept, intercept.cor=intercept.cor,coefficients = coefficients, coefficients.cor=coefficients.cor,
        m.opt = m.opt, m.opt.cor=m.opt.cor,cv.error.matrix = cv.error.matrix, cor.error.matrix=cor.error.matrix,cv.error = cv.error, cor.error=cor.error,
        coefficients.jackknife = coefficients.jackknife))
}

