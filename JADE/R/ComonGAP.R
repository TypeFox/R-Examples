ComonGAP <-
function(A, A.hat)
    {
    # A true mixing matrix
    # A.hat estimated mixing matrix    
    
    if (is.matrix(A) == FALSE) stop("'A' must be a square matrix")
    if (dim(A)[1] != dim(A)[2]) stop("'A' must be a square matrix")
    na.fail(A)
    
    if (is.matrix(A.hat) == FALSE) stop("'A.hat' must be a square matrix")
    if (dim(A.hat)[1] != dim(A.hat)[2]) stop("'A.hat' must be a square matrix")
    na.fail(A.hat)
    
    #A.dash <- A %*% diag(1/sqrt(colSums(A^2)))
    #A.hat.dash <- A.hat%*%diag(1/sqrt(colSums(A.hat^2)))
    
    A.dash <- sweep(A, 2, sqrt(colSums(A^2)), "/")
    A.hat.dash <- sweep(A.hat, 2, sqrt(colSums(A.hat^2)), "/")
    
    D.mat <- solve(A.dash) %*% A.hat.dash
    D.abs.1 <- abs(D.mat)
    D.abs.2 <- D.mat^2
    gap <- sum((rowSums(D.abs.1)-1)^2) + sum((colSums(D.abs.1)-1)^2) + sum(abs(rowSums(D.abs.2)-1)) + sum(abs(colSums(D.abs.2)-1))
    return(gap)
    }
