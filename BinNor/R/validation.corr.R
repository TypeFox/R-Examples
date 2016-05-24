validation.corr <-
function (no.bin, no.nor, prop.vec.bin = NULL, corr.vec = NULL, 
    corr.mat = NULL) 
{
    d = no.bin + no.nor
    validation.bin(no.bin, prop.vec.bin)
    if (is.null(corr.vec) & is.null(corr.mat)) {
        stop("You must specify full correlation matrix OR vector of elements below the diagonal")
    }
    if (!is.null(corr.vec) & !is.null(corr.mat)) {
        corr.mat.from.corr.vec = lower.tri.to.corr.mat(corr.vec, 
            d)
        if (sum(dim(corr.mat.from.corr.vec) == dim(corr.mat)) != 
            2) {
            stop("corr.vec and corr.mat are non-conformable")
        }
        if (sum(corr.mat.from.corr.vec == corr.mat) != (d * d)) {
            stop("Correlation matrix from corr.vec and corr.mat are not the same")
        }
    }
    if (!is.null(corr.vec)) {
        if (length(corr.vec) != (d * (d - 1)/2)) {
            stop("Vector of correlations is misspecified, dimension is wrong!\n")
        }
        if ((min(corr.vec) <= -1) | (max(corr.vec) >= 1)) {
            stop("Correlations must be between -1 and 1!\n")
        }
        corr.mat.from.corr.vec = lower.tri.to.corr.mat(corr.vec, 
            d)
        if (is.positive.definite(corr.mat.from.corr.vec) == FALSE) {
            stop("Specified correlation matrix (from corr.vec) is not positive definite! \n")
        }
        validation.range(no.bin, no.nor, prop.vec.bin, corr.mat.from.corr.vec)
    }
    if (!is.null(corr.mat)) {
        if (dim(corr.mat)[1] != d | dim(corr.mat)[2] != d) {
            stop("Correlation matrix dimension is wrong!\n")
        }
        if (is.positive.definite(corr.mat) == FALSE) {
            stop("Specified correlation matrix is not positive definite! \n")
        }
        if (isSymmetric(corr.mat) == FALSE) {
            stop("Specified correlation matrix is not symmetric! \n")
        }
        validation.range(no.bin, no.nor, prop.vec.bin, corr.mat)
    }
cat("No specification problems are detected for the correlation structure! \n")
}

