matrix.rank <- function (X)
{
    X.sv <- abs(La.svd(X)$d)
    return(sum((X.sv / max(X.sv)) > 1e-09))
}
