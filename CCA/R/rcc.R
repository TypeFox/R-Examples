"rcc" <-
function (X, Y, lambda1, lambda2) 
{
    Xnames <- dimnames(X)[[2]]
    Ynames <- dimnames(Y)[[2]]
    ind.names <- dimnames(X)[[1]]
    Cxx <- var(X, na.rm = TRUE, use = "pairwise") + diag(lambda1, 
        ncol(X))
    Cyy <- var(Y, na.rm = TRUE, use = "pairwise") + diag(lambda2, 
        ncol(Y))
    Cxy <- cov(X, Y, use = "pairwise")
    res <- geigen(Cxy, Cxx, Cyy)
    names(res) <- c("cor", "xcoef", "ycoef")
    scores <- comput(X, Y, res)
    return(list(cor = res$cor, names = list(Xnames = Xnames, 
        Ynames = Ynames, ind.names = ind.names), xcoef = res$xcoef, 
        ycoef = res$ycoef, scores = scores))
}

