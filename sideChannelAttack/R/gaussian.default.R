gaussian.default <-
function (x, y, ...) 
{
    if (!is.matrix(x)) {
        stop("'x' has to be a matrix")
        return(-1)
    }
    if (!is.factor(y)) {
        stop("'y' has to be a factor")
        return(-1)
    }
    if (dim(x)[1] != length(y)) {
        stop("the number of output has to be the same as the number of input")
        return(-1)
    }
    if (length(which(y==0)) > 0) {
        stop("'y' has to contain values strictly greater than 0")
        return(-1)
    }
    X_ = x
    Y_ = y
    CovarianceDuBruitDeLa <- list()
    MoyenneDeLa <- list()
    detCov <- list()
    for (cle in sort(unique(as.numeric(levels(y))[y]))) {
        CovarianceDuBruitDeLa[[cle]] <- matrix(ncol = dim(X_)[2], 
            nrow = dim(X_)[2], data = 0)
        MoyenneDeLa[[cle]] <- matrix(ncol = dim(X_)[2], nrow = 2, 
            data = 0)
        if (length(which(Y_ == cle)) <= 1) {
            stop("the number of data linked to a label has to be strictly greater than 1")
            return(-1)
        }
        CovarianceDuBruitDeLa[[cle]] = cov(as.matrix(X_[which(Y_ == 
            cle), ], ncol = dim(X_)[2]))
        if (det(CovarianceDuBruitDeLa[[cle]]) == 0) {
            CovarianceDuBruitDeLa[[cle]] = matrix(cov.shrink(as.matrix(X_[which(Y_ == 
                cle), ], ncol = dim(X_)[2]), verbose = FALSE), 
                ncol=dim(X_)[2])
        }
	
	if (det(CovarianceDuBruitDeLa[[cle]]) == 0) {
		stop("the determinant of the covariance is zero")
	}

        MoyenneDeLa[[cle]] = apply(matrix(X_[which(Y_ == cle), ], ncol = dim(X_)[2], 
            nrow=length(which(Y_ == cle))),MARGIN = 2, FUN = mean)
        detCov[[cle]] = det(CovarianceDuBruitDeLa[[cle]])
    }
    res = list(mean = MoyenneDeLa, cov = CovarianceDuBruitDeLa, 
        detCov = detCov)
    class(res) <- "gaussian"
    return(res)
}
