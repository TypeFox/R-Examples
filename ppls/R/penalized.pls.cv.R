`penalized.pls.cv` <-
function(X,y,P=NULL,lambda=1,ncomp=NULL,k=5,kernel=FALSE,scale=FALSE)
{
    p<-ncol(X)
    n <- nrow(X)
    if (is.null(ncomp))
        ncomp = min(n - 1, ncol(X))
    lambda = as.vector(lambda)
    if (is.null(P) == TRUE) {
        lambda = 0
    }
    all.folds <- split(sample(1:n), rep(1:k, length = n))
    ntrain = vector(length = k)
    for (i in 1:k) {
        ntrain[i] = n - length(all.folds[[i]])
    }
    ntrain.min = min(ntrain)
    ncomp = min(ncomp, ntrain.min - 1)
    error.cv = matrix(0, length(lambda), ncomp)
    coefficients.jackknife=array(dim=c(p,ncomp,length(lambda),k))
    for (i in seq(k)) {
        omit <- all.folds[[i]]
        Xtrain = X[-omit, , drop = FALSE]
        ytrain = y[-omit]
        Xtest = X[omit, , drop = FALSE]
        ytest = y[omit]
        for (j in 1:length(lambda)) {
            if (is.null(P) == TRUE) {
                Pj = NULL
            }
            if (is.null(P) == FALSE) {
                Pj = lambda[j] * P
            }
            penpls = penalized.pls(Xtrain, ytrain, P = Pj, ncomp,
                kernel = kernel, scale = scale)
                coefficients.jackknife[,,j,i]=penpls$coefficients
            error.cv[j, ] = error.cv[j, ] + length(ytest) * (new.penalized.pls(penpls,
                Xtest, ytest)$mse)
        }
    }
    error.cv = error.cv/n
    value1 = apply(error.cv, 1, min)
    index.lambda=which.min(value1)
    lambda.opt = lambda[index.lambda]
    ncomp.opt = which.min(error.cv[lambda == lambda.opt, ])
    min.ppls = min(value1)
    if (is.null(P) == TRUE) {
        P.opt = NULL
    }
    if (is.null(P) == FALSE) {
        P.opt = lambda.opt * P
    }
    ppls = penalized.pls(X, y, P.opt, ncomp = ncomp.opt, kernel)
    intercept = ppls$intercept[ncomp.opt]
    coefficients = ppls$coefficients[, ncomp.opt]
	outlist=list(error.cv = error.cv,lambda=lambda,ncomp=ncomp,lambda.opt = lambda.opt,index.lambda=index.lambda,
        ncomp.opt = ncomp.opt, min.ppls = min.ppls, intercept = intercept,
        coefficients = coefficients,coefficients.jackknife=coefficients.jackknife)
	class(outlist)="mypls"
    return(outlist)
}