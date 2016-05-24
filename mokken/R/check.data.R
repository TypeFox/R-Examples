"check.data" <-
function(X){
   if (data.class(X) != "matrix" && data.class(X) != "data.frame")
     stop("Data are not matrix or data.frame")
    matrix.X <- as.matrix(X)
    if (any(is.na(matrix.X))) stop("Missing values are not allowed")
    if (mode(matrix.X)!="numeric") stop("Data must be numeric")
    if (any(matrix.X < 0)) stop("All scores should be nonnegative")
    if (any(matrix.X %% 1 !=0)) stop("All scores must be integers")
    matrix.X <- matrix.X - min(matrix.X)
    return(matrix.X)
}
