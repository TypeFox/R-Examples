mrobj <- function (X, m)
{
    return(sum(norme(sweep(X, 2, m))))
}
