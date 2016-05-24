affine.trans <- function(X, A=diag(1, dim(X)[2]), b= rep(0, dim(X)[2]), A.sqrt=FALSE, na.action=na.fail)
    {
    X<-na.action(X)
    
    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric")
    X<-as.matrix(X)
    
    p <- dim(X)[2]
    
    if (!is.logical(A.sqrt)) stop("A.sqrt must be a TRUE or FALSE")
    if (!is.vector(b) | length(b) != p) stop("b must be a vector of length p")
    if (!is.matrix(A) | dim(A)[1] != p | dim(A)[2] !=p ) stop("A must be a p x p matrix")
    
    if (A.sqrt) A <- mat.sqrt(A)
    
    X <- tcrossprod(X,A)
    X <- sweep(X, 2, b, "+")
    
    return(X)
    }
