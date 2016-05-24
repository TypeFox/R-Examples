`delete.na` <-
function (x, alternative = c("less", "greater") )
{
if (alternative == "less")
{

    if (nrow(x) < ncol(x))
        a <- na.omit(x)
    else {
        a <- t(x)
        b <- na.omit(a)
        a <- t(b)
    }
    b <- cbind(a)
    return(b)
}
if (alternative == "greater")
{

    if (nrow(x) > ncol(x))
        a <- na.omit(x)
    else {
        a <- t(x)
        b <- na.omit(a)
        a <- t(b)
    }
    b <- cbind(a)
    return(b)
}
}

