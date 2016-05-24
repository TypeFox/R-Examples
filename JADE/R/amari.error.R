amari.error <-
function(W.hat,A,standardize=F)
    {
    if (is.matrix(W.hat)==FALSE) stop("'W.hat' must be a square matrix")
    if (dim(W.hat)[1]!=dim(W.hat)[2]) stop("'W.hat' must be a square matrix")
    na.fail(W.hat)
    
    if (is.matrix(A)==FALSE) stop("'A' must be a square matrix")
    if (dim(A)[1]!=dim(A)[2]) stop("'A' must be a square matrix")
    na.fail(A)
    
    if (dim(A)[1]!=dim(W.hat)[1]) stop("'W.hat' and 'A' must have the same dimensions")
    if (dim(A)[1]<2) stop("'W.hat' and 'A' must be at least 2x2 matrices")
    
    if (!is.logical(standardize)) stop("'standardize' must be logical")
    
    if (standardize==TRUE)
        {
        W.hat <- .standard.B(W.hat)
        A <- solve(.standard.B(solve(A)))
        }
    
    P.abs <- abs(W.hat %*% A)
    k <- dim(W.hat)[1]
    row.max <- apply(P.abs,1,max)
    col.max <- apply(P.abs,2,max)
    E.1 <- sweep(P.abs,1,row.max,"/")
    E.2 <- sweep(P.abs,2,col.max,"/")
    A.error <- 1/(2*k*(k-1)) * (sum(rowSums(E.1)-1) + sum(colSums(E.2)-1))
    return(A.error)
    }
