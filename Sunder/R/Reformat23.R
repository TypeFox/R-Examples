Reformat23 <-
function (x) 
{
    if (length(dim(x)) != 2) 
        stop("x must be 2d")
    if (class(x) == "data.frame") 
        x = data.matrix(x)
    nind = nrow(x)
    nc = ncol(x)
    nloc = nc/2
    y = tapply(as.vector(x), rep(1:nloc, each = 2 * nind), function(x) as.numeric(as.factor(x)))
    ym = matrix(unlist(y), nind)
    nal = as.vector(sapply(y, max, na.rm = TRUE))
    nalM = max(nal)
    g = array(dim = c(nind, nloc, nalM), 0)
    c1 = as.vector(row(ym))
    c2 = ceiling(as.vector(col(ym))/2)
    c3 = as.vector(ym)
    qq = table(nloc * nind * (c3 - 1) + nind * (c2 - 1) + c1)
    g[as.numeric(names(qq))] = qq
    return(g)
}
