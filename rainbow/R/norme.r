norme <- function (X)
{
    return(sqrt(rowSums(X^2, na.rm = TRUE)))
}
