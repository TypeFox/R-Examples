rinvchisq <- function (n, df, scale = 1/df) 
{
    if ((length(scale) != 1) & (length(scale) != n)) 
        stop("scale should have the same length as n")
    if (df <= 0) 
        stop("df must be greater than zero")
    if (any(scale <= 0)) 
        stop("scale must be greater than zero")
    (df * scale)/rchisq(n, df = df)
}