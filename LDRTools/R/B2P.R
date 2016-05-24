B2P <-
function(x)
    {
    x <- as.matrix(x)
    x %*% tcrossprod(solve(crossprod(x)), x) 
    }
