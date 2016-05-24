root.matrix <- function(X)
{
    if((ncol(X)==1)&&(nrow(X)==1)) return(sqrt(X))
    else
    {
        X.eigen <- eigen(X, symmetric=TRUE)
        if(any(X.eigen$values < 0)) stop("matrix is not positive semidefinite")
        sqomega <- sqrt(diag(X.eigen$values))
        V <- X.eigen$vectors
        return(V %*% sqomega %*% t(V))
    }
}

solveCrossprod <- function(X, method = c("qr", "chol", "solve")) {
  switch(match.arg(method),
    "qr" = chol2inv(qr.R(qr(X))),
    "chol" = chol2inv(chol(crossprod(X))),
    "solve" = solve(crossprod(X)))
}


