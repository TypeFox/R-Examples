validation.corr <-
function(n.P, n.B, n.O, corr.vec = NULL, corr.mat = NULL) {
        
    if (missing(n.P) == TRUE) {
       stop("Number of Poisson variables is not specified!")
    }  else
    if (missing(n.B) == TRUE) {
       stop("Number of binary variables is not specified!")
    }  else
    if (missing(n.O) == TRUE) {
       stop("Number of ordinal variables is not specified!")
    }  else
    if (!missing(n.P) && !missing(n.B) && !missing(n.O)) {
       
        if ((n.P < 0) | (floor(n.P) != n.P))    {
        stop("Number of Poisson variables must be a non-negative integer!")
        } else 
        if ((n.B < 0) | (floor(n.B) != n.B))    {
        stop("Number of binary variables must be a non-negative integer!")
        } else
        if ((n.O < 0) | (floor(n.O) != n.O))    {
        stop("Number of ordinal variables must be a non-negative integer!")
        } else

        d = n.P + n.B + n.O
    } #if


    if (is.null(corr.mat) & is.null(corr.vec))   {
        stop("You must specify full correlation matrix OR vector of elements below the diagonal!")
    } #if

    if (!is.null(corr.mat) & !is.null(corr.vec)) {
        corr.mat.from.corr.vec=diag(1,d)
        corr.mat.from.corr.vec[lower.tri(corr.mat.from.corr.vec)]=corr.vec
        corr.mat.from.corr.vec=corr.mat.from.corr.vec+t(corr.mat.from.corr.vec)-diag(1,d)
        if (sum(dim(corr.mat.from.corr.vec) == dim(corr.mat)[1]) !=2) {
        stop("corr.vec and corr.mat are non-conformable!")
        }#if
        if (sum(corr.mat.from.corr.vec == corr.mat) != (d * d))   {
            stop("Correlation matrix from corr.vec and corr.mat are not the same!")
        }#if
    } #if

    if (!is.null(corr.vec)) {
        if (length(corr.vec) != (d * (d - 1)/2)) {
            stop("Vector of correlations is misspecified, dimension is wrong!\n")
        } #if
        if ((min(corr.vec) <= -1) | (max(corr.vec) >= 1)) {
            stop("Correlations must be between -1 and 1!\n")
        } #if
        corr.mat.from.corr.vec=diag(1,d)
        corr.mat.from.corr.vec[lower.tri(corr.mat.from.corr.vec)]=corr.vec
        corr.mat.from.corr.vec=corr.mat.from.corr.vec+t(corr.mat.from.corr.vec)-diag(1,d)
        if (is.positive.definite(corr.mat.from.corr.vec) == FALSE) {
            stop("Specified correlation matrix (from corr.vec) is not positive definite! \n")
        } #if
      
    }#if

    if (!is.null(corr.mat)) {
        if (dim(corr.mat)[1] != d | dim(corr.mat)[2] != d) {
            stop("Correlation matrix dimension is wrong!\n")
        }#if
        if (is.positive.definite(corr.mat) == FALSE) {
            stop("Specified correlation matrix is not positive definite! \n")
        }#if
        if (isSymmetric(corr.mat) == FALSE) {
            stop("Specified correlation matrix is not symmetric! \n")
        }#if
        
    }#if
return(TRUE)
}
