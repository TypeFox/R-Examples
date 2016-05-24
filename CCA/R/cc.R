"cc" =
function (X, Y) 
{
    Xnames = dimnames(X)[[2]]
    Ynames = dimnames(Y)[[2]]
    ind.names = dimnames(X)[[1]]
    res = rcc(X, Y, 0, 0)
    return(res)
}

