"comput" =
function (X, Y, res) 
{
    X.aux = scale(X, center=TRUE, scale=FALSE)
    Y.aux = scale(Y, center=TRUE, scale=FALSE)
    X.aux[is.na(X.aux)] = 0
    Y.aux[is.na(Y.aux)] = 0

    xscores = X.aux%*%res$xcoef
    yscores = Y.aux%*%res$ycoef

    corr.X.xscores = cor(X, xscores, use = "pairwise")
    corr.Y.xscores = cor(Y, xscores, use = "pairwise")
    corr.X.yscores = cor(X, yscores, use = "pairwise")
    corr.Y.yscores = cor(Y, yscores, use = "pairwise")

    return(list(xscores = xscores, yscores = yscores, 
        corr.X.xscores = corr.X.xscores, 
        corr.Y.xscores = corr.Y.xscores, 
        corr.X.yscores = corr.X.yscores, 
        corr.Y.yscores = corr.Y.yscores))
}

