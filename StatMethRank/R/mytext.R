mytext <- function (yval, digits, n, use.n)
{
    if (use.n)
        return(paste0(yval, "\nn=", n))
    else
        return(paste0(yval))
}