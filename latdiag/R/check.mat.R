check.mat <- function (X)   # modified from mokken
{
    if (data.class(X) != "matrix" && data.class(X) != "data.frame") 
        stop("mat should be matrix or data.frame")
    matrix.X <- as.matrix(X)
    if (any(is.na(matrix.X))) 
        stop("Missing values not allowed in mat")
    if (mode(matrix.X) != "numeric") 
        stop("mat must be numeric")
    if (any(matrix.X < 0)) 
        stop("All scores should be nonnegative")
    if (any(matrix.X %% 1 != 0)) 
        stop("All scores must be integers")
    matrix.X <- matrix.X - min(matrix.X)
   if(length(table(matrix.X)) != 2)
      stop("mat must be binary")
   matrix.X <- matrix.X / max(matrix.X)
    return(matrix.X)
}
