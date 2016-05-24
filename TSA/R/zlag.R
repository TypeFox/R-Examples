`zlag` <-
function (x, d = 1) 
{
    if (d != as.integer(d) || d < 0) 
        stop("d must be a non-negative integer")
    if (d == 0) 
        return(x)
    else return(c(rep(NA, d), rev(rev(x)[-(1:d)])))
}
