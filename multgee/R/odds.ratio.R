odds.ratio <-
function (x) 
{
    dimx <- nrow(x)
    x[-1, -1] * x[-dimx, -dimx]/x[-dimx, -1]/x[-1, -dimx]
}

